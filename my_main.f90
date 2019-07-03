program my_main

  use HDF5
  use MPI

  implicit none

  !include 'mpi.h'

  ! Numbers of things:
  integer :: nvalence  
  integer :: nconduction 
  integer :: nbands
  integer :: nbands_to_write
  integer :: ncols = 5 ! in data file
  integer :: npools 
  integer :: nvownmax 
  integer :: ncownmax

  ! Distribution subroutine stuff:
  integer :: mypool, mypoolrank, myipe
  integer :: nvownactual, ncownactual
  integer, allocatable  :: global_pairowner(:,:)
  logical, allocatable  :: doiownv(:)
  logical, allocatable  :: doiownc(:)
  logical, allocatable  :: does_it_ownv(:,:)
  logical, allocatable  :: does_it_ownc(:,:)
  integer, allocatable  :: global_nvown(:)
  integer, allocatable  :: global_ncown(:)
  integer, allocatable  :: indexv(:)
  integer, allocatable  :: global_indexv(:,:)
  integer, allocatable  :: indexc(:)
  integer, allocatable  :: invindexv(:)
  integer, allocatable  :: invindexc(:)
  integer :: i, j, k, ic, iv, ipool, ipe

  ! Temporary variables for distribution:
  integer, allocatable :: global_pairowner_temp(:,:)
  integer, allocatable :: global_nvown_temp(:)
  integer, allocatable :: global_ncown_temp(:)
  integer, allocatable :: global_indexv_temp(:,:)

  ! MPI stuff:
  integer :: mpierror
  integer :: comm, info
  integer :: inode, npes, an_id
  integer :: file_write, ierr!, status(MPI_STATUS_SIZE), sender
  integer, parameter :: return_data_tag = 2001


  ! Inclusion array stuff:
  integer, dimension(4,2) :: incl_array
  integer, allocatable :: incl_array_v(:,:)
  integer, allocatable :: incl_array_c(:,:)
  integer, allocatable :: my_incl_array_v(:,:), my_incl_array_c(:,:)

  ! HDF5 stuff:
  character(len=10), parameter :: filename = "sds_row.h5"
  character(len=8), parameter :: dsetname = "IntArray"
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: filespace
  integer(HID_T) :: memspace
  integer(HID_T) :: dataspace
  integer(HID_T) :: plist_id
  integer(HSIZE_T), dimension(2) :: dimsf
  integer(HSIZE_T), dimension(2) :: dimsm
  integer(HSIZE_T), dimension(2) :: count, count_out
  integer(HSIZE_T), dimension(2) :: offset, offset_out
  integer(HSIZE_T), dimension(2) :: stride
  integer(HSIZE_T), dimension(2) :: block
  integer, allocatable :: data(:,:) ! Data to write
  integer, allocatable :: my_bands_v(:,:) ! Data to read back in
  integer, allocatable :: my_bands_c(:,:) ! Data to read back in
  integer :: error ! hdf5 error flag

  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  call MPI_INIT(mpierror)
  call MPI_COMM_SIZE(comm, npes, mpierror)
  call MPI_COMM_RANK(comm, inode, mpierror)

  call h5open_f(error)

!   Define a dummy inclusion array
  incl_array(1,1) = 2
  incl_array(1,2) = 2
  incl_array(2,1) = 4
  incl_array(2,2) = 6
  incl_array(3,1) = 8
  incl_array(3,2) = 9
  incl_array(4,1) = 12
  incl_array(4,2) = 14
  nvalence = 5
  nconduction = 4
  nbands = nvalence + nconduction
  nbands_to_write = incl_array(4,2) + 1
  ! rows included: 2, 4, 5, 6, 8, 9, 12, 13, 14
  ! rows_v incl: 2, 4, 5, 6, 8
  ! rows_c incl: 9, 12, 13, 14

  ! Make the dummy datafile to read
  call make_file()
  call MPI_Barrier(comm, mpierror) ! wait for task 0 to finish
  ! Split the inclusion_array into valence and conduction parts
  call make_vc_incl_array(incl_array, nvalence, nconduction, &
    incl_array_v, incl_array_c)

  if (inode .eq. 0) then
    write(*,*) "incl_array: "
    do i = 1, size(incl_array,1)
      write(*,*) (incl_array(i, j), j =1,2)
    end do
    write(*,*) "incl_array_v: "
    do i = 1, size(incl_array_v,1)
      write(*,*) (incl_array_v(i, j), j =1,2)
    end do
    write(*,*) "incl_array_c: "
    do i = 1, size(incl_array_c,1)
      write(*,*) (incl_array_c(i, j), j =1,2)
    end do
  end if
  ! Distribute bands evenly to the mpi tasks
  call my_distribution()

  ! Make the valence and conduction local inclusion arrays
  call make_my_incl_array(incl_array_v, doiownv, my_incl_array_v, nvownactual)
  call make_my_incl_array(incl_array_c, doiownc, my_incl_array_c, ncownactual)
!  call read_bands

   do k = 1, npes
     if (inode .eq. k-1) then
       write(*,*) "task:  ", inode, "my_incl_array_v: "
       do i = 1, nvownactual
         write(*,*) (my_incl_array_v(i,j), j = 1, 2)
       end do
       write(*,*) "task:  ", inode, "my_incl_array_c: "
       do i = 1, ncownactual
         write(*,*) (my_incl_array_c(i,j), j = 1, 2)
       end do
     end if
     call MPI_Barrier(comm, mpierror)
   end do

  allocate(my_bands_v(nvownactual, ncols))
  allocate(my_bands_c(ncownactual, ncols))

  call read_bands(my_incl_array_v, nvownactual, ncols, my_bands_v)
  call read_bands(my_incl_array_c, ncownactual, ncols, my_bands_c)

  do k = 1, npes
    if (inode .eq. k-1) then
      write(*,*) "task:  ", inode, " my_bands_v: "
      do i = 1, nvownactual
        write(*,*) (my_bands_v(i,j), j = 1, ncols)
      end do
      write(*,*) "task:  ", inode, " my_bands_c: "
      do i = 1, ncownactual
        write(*,*) (my_bands_c(i,j), j = 1, ncols)
      end do
    end if
    call MPI_Barrier(comm, mpierror)
  end do

  deallocate(global_pairowner)
  deallocate(doiownv)
  deallocate(doiownc)
  deallocate(does_it_ownv)
  deallocate(does_it_ownc)
  deallocate(global_nvown)
  deallocate(global_ncown)
  deallocate(indexv)
  deallocate(global_indexv)
  deallocate(indexc)
  deallocate(invindexv)
  deallocate(invindexc)

  deallocate(incl_array_v)
  deallocate(incl_array_c)
  !deallocate(my_incl_array_v)
  deallocate(my_incl_array_c)

  !deallocate(data)
  !deallocate(data_out)

  call h5close_f(error)
  call MPI_FINALIZE(mpierror)

  contains

    subroutine my_distribution()

      call createpools(nvalence, nconduction, npes, npools, nvownmax, ncownmax)

     ! if (npools .gt. npes) then
        if (mod(nvalence, npools) .eq. 0) then
          nvownmax = nvalence / npools
        else
          nvownmax = (nvalence/npools) + 1
        end if
    
        if (mod(nconduction, (npes/npools)) .eq. 0) then
          ncownmax = (nconduction) / (npes/npools)
        else
          ncownmax = (nconduction) / (npes/npools) + 1
        end if
     ! end if

      ! Make allocations:
      allocate( global_pairowner(nvalence, nconduction))
      global_pairowner(:,:) = 0
      allocate( doiownv(nvalence))
      doiownv(:) = .false.
      allocate( doiownc(nconduction))
      doiownc(:) = .false.
      allocate( does_it_ownv(nvalence, npes))
      does_it_ownv(:,:) = .false.
      allocate( does_it_ownc(nconduction, npes))
      does_it_ownc(:,:) = .false.
      allocate( global_nvown(npes))
      global_nvown(:) = 0
      allocate( global_ncown(npes))
      global_ncown(:) = 0
      allocate( indexv(nvalence))
      indexv(:) = 0
      allocate( global_indexv(nvalence, npes))
      global_indexv(:,:) = 0
      allocate( indexc( nconduction))
      indexc(:) = 0
      allocate( invindexv( nvownmax))
      invindexv(:) = 0
      allocate( invindexc( ncownmax))
      invindexc(:) = 0
      ! try to do this after lunch
    
      nvownactual = 0
      ncownactual = 0
    
      mypool = inode/(npes/npools)
      mypoolrank = mod(inode, (npes/npools))
      myipe = inode + 1
!      write(*,*) "My ID: ", inode, "My pool: ", mypool, "My pool rank: ", mypoolrank
    
      do iv = 1, nvalence
        ipool = (iv - 1)/nvownmax
    
        if (mypool .eq. ipool .and. inode .lt. npes) then
    
          nvownactual = nvownactual+1
          global_nvown(myipe) = nvownactual
          indexv(iv) = nvownactual
          global_indexv(iv, myipe) = indexv(iv)
          invindexv(nvownactual) = iv
          doiownv(iv) = .true.
    
          do ic = 1, nconduction
            if ( (ic - 1)/ ncownmax == mypoolrank) then
              if (nvownactual .eq. 1) then
                 ncownactual = ncownactual + 1
                 global_ncown(myipe) = ncownactual
                 invindexc(ncownactual) = ic
                 indexc(ic) = ncownactual
                 doiownc(ic) = .true.
               end if
               global_pairowner(iv, ic) = myipe
             end if
         end do ! end do ic = 1, nconduction
    
        endif ! end if my pool the current pool
      end do ! end do iv = 1, nvalence
    
      if (inode .lt. npes) then
        does_it_ownv(:, inode+1) = doiownv(:)
        does_it_ownc(:, inode+1) = doiownc(:)
      end if ! inode < npes
    
    
      do ipe=0, npes-1
        call MPI_Bcast(does_it_ownv(1, ipe+1), nvalence, &
          MPI_LOGICAL, ipe, MPI_COMM_WORLD, mpierror)
        call MPI_Bcast(does_it_ownc(1, ipe+1), nconduction, &
          MPI_LOGICAL, ipe, MPI_COMM_WORLD, mpierror)
      end do !ipe = 0, npes-1
    
      allocate(global_pairowner_temp(nvalence, nconduction))
      call MPI_ALLREDUCE(global_pairowner(1,1), global_pairowner_temp(1,1), &
        nvalence*nconduction, MPI_INTEGER, MPI_SUM, &
        MPI_COMM_WORLD, mpierror)
      global_pairowner = global_pairowner_temp
      deallocate(global_pairowner_temp)
    
      allocate(global_nvown_temp(npes))
      call MPI_ALLREDUCE(global_nvown, global_nvown_temp, npes, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_WORLD, mpierror)
      global_nvown = global_nvown_temp
      deallocate(global_nvown_temp)
    
      allocate(global_ncown_temp(npes))
      call MPI_ALLREDUCE(global_ncown, global_ncown_temp, npes, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_WORLD, mpierror)
      global_ncown = global_ncown_temp
      deallocate(global_ncown_temp)
    
      allocate(global_indexv_temp(nvalence, npes))
      call MPI_ALLREDUCE(global_indexv(1,1), global_indexv_temp(1,1), &
        npes*nvalence, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
      global_indexv = global_indexv_temp
      deallocate(global_indexv_temp)
    
      call MPI_BARRIER(comm, mpierror)
    
      if ( inode .eq. 0) then

        write(*,*) "Does_it_ownv:"
        do i = 1, nvalence
          write(*,*) (does_it_ownv(i, j), j = 1, npes)
        end do
        
        write(*,*) "Does_it_ownc:"
        do i = 1, nconduction
          write(*,*) (does_it_ownc(i, j), j = 1, npes)
        end do
        
      end if ! if inode == 0

    end subroutine my_distribution

      
    subroutine createpools(npq, nsq, npess, npoolsout, npqownmaxout, nsqownmaxout)

      integer, intent(in) :: npq, nsq, npess
      integer, intent(out) :: npoolsout, npqownmaxout, nsqownmaxout
      integer :: nmemmin, npoolss, nsqownmax, npqownmax, nmem, npes_per_npoolss

      nmemmin = npq + nsq + 1
      npoolsout = 0
      npqownmaxout = 0
      nsqownmaxout = 0

      do npoolss = 1, min(npq, npess)
        npqownmax = (npq + npoolss - 1) / npoolss
        npes_per_npoolss = npess/npoolss
        nsqownmax = (nsq + npes_per_npoolss - 1) / npes_per_npoolss

        nmem = npqownmax + nsqownmax

        if (nmem .lt. nmemmin) then
          nmemmin = nmem
          npoolsout = npoolss
          npqownmaxout = npqownmax
          nsqownmaxout = nsqownmax
        end if

      end do

      return

    end subroutine createpools

    subroutine make_my_incl_array(incl_array, do_i_own, my_incl_array,&
      nownactual)

      integer, intent(in) :: incl_array(:,:)
      integer, intent(in) :: nownactual
      integer, intent(out), allocatable :: my_incl_array(:, :)
      logical, intent(in) :: do_i_own(:)
    
      integer :: current_last_incl_row
      integer :: nrows, ncols
      !integer :: my_count = 0 ! for allocating my_incl_array
      integer :: irow ! row indexer of my_incl_array
      integer :: next_row
      integer :: i_ia  ! "(i)ndex of (i)nclusion (a)rray" 
      integer :: i, j, k
      integer :: init 
      
      nrows = SIZE(incl_array, 1) ! first dim of incl_array
      ncols = 2 ! Fixed by def of inclusion array
    
      ! the task does not own any bands, then an error will occur
      ! When trying to read the actual bands in
      ! but by setting to -1, we ensure that it gets skipped over
      ! in the reading routine
      allocate(my_incl_array(nownactual, ncols)) ! allocate worst case scenario
      my_incl_array = -1 ! -1 if row is extra
      init = 0
    
      next_row = incl_array(1,1)
      i_ia = 1
      irow = 2
      do i = 1, size(do_i_own)
        if (do_i_own(i)) then
          if ( init .eq. 0) then
            ! if I am the first addition to my_incl_array, then just add me with
            ! no further logic:
            my_incl_array(1, 1) = next_row
            my_incl_array(1, 2) = next_row
            init = init + 1
          else
            ! Must include i in my_incl_array, but have to figure out if i
            ! needs its own row, or should be added to the end of the current row
            ! grab the bottom right value in my_incl_array
            current_last_incl_row = my_incl_array(irow - 1, 2)
            ! if the next value is adjacent, then just increment the bottom right
            ! value to the next row.
            if ( next_row .eq. current_last_incl_row + 1 ) then
              my_incl_array(irow - 1, 2) = my_incl_array(irow - 1, 2) + 1
  !            my_incl_array(irow, 2) = my_incl_array(irow, 2) + 1
            else ! if not adjacent, then the band is not in the current cluster, so
              ! start a new one 
              irow = irow + 1
              my_incl_array(irow - 1, 1) = next_row
              my_incl_array(irow - 1, 2) = next_row
            end if ! (next_row...)
          end if ! init .eq. 0
        end if ! do_i_own(i)
        if ( i .ne. size(do_i_own)) then
          if (incl_array(i_ia, 2) .gt. next_row) then
            next_row = next_row + 1
          else
            next_row = incl_array(i_ia+1, 1)
           i_ia = i_ia + 1
          end if ! incl_array(i_ia) ... 
        end if ! i .ne. size(do_i_own)
      end do
    
    ! write(*,*) " "
    !  do i = 1, nownactual
    !    write(*,*) (inode, my_incl_array(i, j), j = 1, ncols)
    !  end do
    !  write(*,*) " "
    
      end subroutine make_my_incl_array

      subroutine make_vc_incl_array(incl_array, nvalence, nconduction, &
          incl_array_v, incl_array_c)

        integer, intent(in) :: incl_array(:,:)
        integer, intent(in) :: nvalence
        integer, intent(in) :: nconduction    
        integer, allocatable, intent(out) :: incl_array_v(:,:)
        integer, allocatable, intent(out) :: incl_array_c(:,:)

        integer :: vcount = 0
        integer :: j = 0
        integer :: nrows
        integer :: find_v
        integer :: crows
        integer :: k

        nrows = size(incl_array, 1)
        
        do while (vcount .lt. nvalence)
          j = j + 1
          vcount = vcount + incl_array(j, 2) - incl_array(j, 1) + 1
        end do
        
        crows = nrows - j
        find_v = incl_array(j, 2)
        
        if (vcount .eq. nvalence) then
          allocate(incl_array_c(crows, 2))
          allocate(incl_array_v(j, 2))
          incl_array_v = incl_array(1:j, :)
          incl_array_c = incl_array(j+1:nrows, :)
        else 
          do while (vcount .ne. nvalence)
            vcount = vcount - 1
            find_v = find_v - 1
          end do
          allocate(incl_array_v(j,2))
          allocate(incl_array_c(crows+1 , 2))
        
          incl_array_v = incl_array(1:j, :)
          incl_array_c = incl_array(j:nrows, :)
          incl_array_v(j, 2) = find_v 
          incl_array_c(1, 1) = find_v + 1
        end if
        
     !   write(*,*) "incl_array_v: "
     !   do i = 1, j
     !     write(*,*) (incl_array_v(i, k), k = 1, 2)
     !   end do
     !   
     !   write(*,*) "incl_array_c: "
     !   write(*,*)
     !   do i = 1, size(incl_array_c, 1)
     !     write(*,*) (incl_array_c(i, k), k = 1,2)
     !   end do

      end subroutine make_vc_incl_array
     

      subroutine make_file()

        allocate(data( nbands_to_write, ncols) )

        if (inode .eq. 0) then
          do i = 1, nbands_to_write
            do j = 1, ncols
              data(i, j) = i
            end do
          end do

          dimsf(1) = nbands_to_write
          dimsf(2) = ncols

          call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
          call h5screate_simple_f(2, dimsf, dataspace, error)
          call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
            dset_id, error)
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsf, error)
          call h5sclose_f(dataspace, error)
          call h5dclose_f(dset_id, error)
          call h5fclose_f(file_id, error)

          do i = 1, nbands_to_write
            write(*,*) (data(i,j), j = 1, ncols)
          end do

          file_write = 1

        end if

      end subroutine make_file

      subroutine read_bands(my_incl_array, row_dim, col_dim, data_out) 
        integer, allocatable, intent(in) :: my_incl_array(:,:) ! v or c
        integer, intent(in) :: row_dim
        integer, intent(in) :: col_dim
        integer, intent(inout) :: data_out(:,:)

        integer :: incl_array_nrows

        incl_array_nrows = size(my_incl_array, 1)

        dimsm(1) = row_dim ! ncownactual or nvownactual
        dimsm(2) = col_dim

        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, comm, info, mpierror)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, & 
          error, access_prp=plist_id)
        call h5pclose_f(plist_id, error)

        call h5dopen_f(file_id, dsetname, dset_id, error)
        call h5dget_space_f(dset_id, dataspace, error)
        call h5screate_simple_f(2, dimsm, memspace, error)

        offset(2) = 0
        offset_out(2) = 0
        count(2) = ncols
        count_out(2) = ncols

        if (row_dim .ne. 0) then

          offset(1) = my_incl_array(1,1) - 1
          count(1) = my_incl_array(1,2) - my_incl_array(1,1) + 1
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, &
            count, error)
          
          do i = 2, incl_array_nrows
            if (my_incl_array(i, 1) .ne. -1) then
              offset(1) = my_incl_array(i,1) - 1
              count(1) = my_incl_array(i,2) - my_incl_array(i, 1) + 1
              if (inode .eq. 0) then
              end if
              call H5sselect_hyperslab_f(dataspace, H5S_SELECT_OR_F, offset, &
                count, error)
            end if
          end do

          offset_out(1) = 0
          count_out(1) = my_incl_array(1,2) - my_incl_array(1,1) + 1

          call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_out, &
            count_out, error)
          
          offset_out(1) = count_out(1) 
         ! offset_out(1) = 0
          do i = 2, incl_array_nrows
            if (my_incl_array(i, 1) .ne. -1) then
              count_out(1) = my_incl_array(i, 2) - my_incl_array(i, 1) + 1
              if (inode .eq. 0) then
              end if
              call h5sselect_hyperslab_f(memspace, H5S_SELECT_OR_F, offset_out, &
                count_out, error)
              offset_out(1) = offset_out(1) + count_out(1)
            end if
          end do

        else
          count(1) = 0
          offset(1) = 0
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, &
            count, error)
          call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, &
            count, error)
        end if ! row dim ne 0
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        call H5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, dimsm, error, &
          memspace, dataspace, xfer_prp=plist_id)

        call h5pclose_f(plist_id, error)
        call h5sclose_f(memspace, error)
        call h5sclose_f(dataspace, error)
        call h5dclose_f(dset_id, error)
        call h5fclose_f(file_id, error)

      end subroutine read_bands

end program my_main
