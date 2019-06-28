/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* We do support the IEC 559 math functionality, real and complex.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
!>=========================================================================
!!
!! Module:
!!
!! (1) wfn_io_hdf5_m Originally by JIM Last Modified 4/25/2012 (JIM)
!!
!! Routines to read and write wavefunctions in 1 format.
!! The code is generated through repeated inclusion of a file with
!! different preprocessor definitions each time. You are not expected to
!! understand this. Consult the resulting .p.f file for clarity.
!!
!!=========================================================================
!The following macro puts any point/array in the [-0.5, 0.5) range:
!The following macro puts any point/array in the [0, 1) range:
!Integer division of a/b rounded up*/
!Rounds a up to the smallest multiple of b*/
! disable Fortran 1 pragmas if not -DOMP*/
! note: C standard does not permit $ in identifiers, however this seems acceptable
! as an extension, for all versions of cpp I tried. --DAS
! truncate spaces in string
!#!define TRUNC(s) trim(adjustl(s))
! Sun compiler has a length limit of 132 characters and won`t support these macros
! No checking for faster performance, if not in debug mode
! Use this instead of the intrinsic 'deallocate' for pointers
! Use this instead of the intrinsic 'deallocate' for arrays
!the TOSTRING macro converts a macro into a string
! deprecated identifiers
! Created Sept 2011 by DAS.
! Define characteristics of various compilers, via compiler symbols (e.g. -DGNU)
! to be used directly from the arch.mk files, and then defining what we need to do
! for that compiler via the symbols for various properties (e.g. NOSIZEOF).
! Ideally, to support a new compiler, one need only change this file, adding a
! new block to define what -DNEWCOMPILER would mean.
! NOTE: of course, Makefile-level issues still need to be handled in common-rules.mk
! very ancient version may require NOSIZEOF
! open64 is very similar to path, it is an open-sourced version of it
! omp_lib.f90 needed to do OpenMP, see common-rules.mk.
! path before 4.0.9 lacks SIZEOF
! both open64 and path die on fseek with:
!lib-5002 : UNRECOVERABLE library error
! This FFIO request is not supported.
!
!Encountered during a GETPOS on unit 8
! cce 7.4.4 and before support sizeof for intrinsic types, but need NOSIZEOF_TYPE
! cce 8.0.0 and later do not allow sizeof for multidimensional arrays, requiring us
! to turn sizeof off everywhere. Why would Cray do this?
! It is considered a bug in OPEN64 that sizeof will not work in our code.
! Fortran 2003 prefers this as a statement, not an intrinsic
! on some platforms there is a different return value for sizeof if build is 64-bit
! name of routine to get name of host program is running on
! module required for hostnm routine to be usable
! note: interfaces are split into two lines because ifort for some reason decrees of
! 'end interface' : "END statement must be only statement on line".
! HOSTNAMEINT enables interface for hostnm routine in intrinsics_m
! how to get the cpu time in seconds
! interface required for mclock routine (timing) to be usable
! interface required for iargc routine (command-line arguments) to be usable
! ftell gives you the current location in a file, to fseek back to it
! if no fseek, ftell is useless
! #warning This compiler does not support fseek.
! fseek returns to a location in a file, bookmarked by ftell. G95 lacks it
! intrinsic module for OpenMP. external for Open64 (see common-rules.mk). NAG and G95 do not support OpenMP
! using a global var here to avoid need for conditional local declaration
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
module wfn_io_hdf5_m
  use global_m
  use hdf5_io_m
  use hdf5
  implicit none
  private
  public :: &
    read_hdf5_header_type , &
    write_hdf5_header_type , &
    read_hdf5_gvectors , &
    write_hdf5_gvectors , &
    read_hdf5_wfn_gvectors , &
    write_hdf5_wfn_gvectors , &
    read_hdf5_band_real , &
    write_hdf5_band_real , &
    read_hdf5_band_complex , &
    write_hdf5_band_complex , &
    read_hdf5_band , &
    write_hdf5_band , &
    read_hdf5_bands_block , &
    setup_hdf5_mf_file , &
    setup_hdf5_wfn_file , &
    read_hdf5_mf_header , &
    write_hdf5_mf_header
  interface read_hdf5_band
    module procedure read_hdf5_band_real, read_hdf5_band_complex
  end interface
  interface write_hdf5_band
    module procedure write_hdf5_band_real, write_hdf5_band_complex
  end interface
contains
!===============================================================================
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
subroutine read_hdf5_header_type(sFileName, sheader, iflavor, kp, gvec, syms, crys)
  character(len=*), intent(in) :: sFileName
  character(len=3), intent(inout) :: sheader
  !> define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(inout) :: iflavor
  type(kpoints), intent(out) :: kp
  type(gspace), intent(out) :: gvec
  type(symmetry), intent(out) :: syms
  type(crystal), intent(out) :: crys
  ! set values based on epsilon calculation
  logical :: is_get=.false.
  logical :: wfnflag=.true.
 
  if (peinf%inode == 0) then
    call read_info(TRUNC(sFileName),iflavor)
    call read_kpoints(TRUNC(sFileName),kp)
    call read_gspace(TRUNC(sFileName),gvec)
    call read_symmetry(TRUNC(sFileName),syms)
    call read_crystal(TRUNC(sFileName),crys)
  endif
  if (peinf%npes > 1) then
    call MPI_BCAST(iflavor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(is_get, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kp%nspin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(kp%nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(gvec%ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(gvec%ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(gvec%FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(syms%tnp, 3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    if (wfnflag) then
      call MPI_BCAST(kp%nrk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%mnband, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ngkmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%kgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%shift, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    endif
  endif
  if(peinf%inode > 0 ) then
    allocate(crys%atyp (crys%nat))
    allocate(crys%apos (3, crys%nat))
    if (wfnflag) then
      allocate(kp%ngk (kp%nrk))
      allocate(kp%w (kp%nrk))
      allocate(kp%rk (3, kp%nrk))
      allocate(kp%ifmin (kp%nrk, kp%nspin))
      allocate(kp%ifmax (kp%nrk, kp%nspin))
      allocate(kp%el (kp%mnband, kp%nrk, kp%nspin))
      allocate(kp%occ (kp%mnband, kp%nrk, kp%nspin))
    endif
  endif
  if (peinf%npes > 1) then
    call MPI_BCAST(crys%atyp, crys%nat, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(crys%apos, 3*crys%nat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    if (wfnflag) then
      call MPI_BCAST(kp%ngk, kp%nrk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%w, kp%nrk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%rk(1,1), 3*kp%nrk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ifmin(1,1), size(kp%ifmin), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%ifmax(1,1), size(kp%ifmax), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%el(1,1,1), size(kp%el), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
      call MPI_BCAST(kp%occ(1,1,1), size(kp%occ), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
    endif
  endif
 
end subroutine read_hdf5_header_type
subroutine read_info(sFileName, iflavor)
  character(len=*), intent(in) :: sFileName
  !> define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(inout) :: iflavor
  integer(HID_T) :: hidFile
  integer :: iError, iflavor2
  character(len=16) :: sflavor
  character(len=128) :: routine_name
  logical :: exists
 
  call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
  iflavor2 = iflavor
  if (iflavor==0) iflavor2 = 1
  routine_name = "read_info"
  ! FHJ: Reading flavor: need some special logic b/c of iflavor==-1
  if (iflavor<-1 .or. iflavor>2) then
    write(sflavor,'(i0)') iflavor
    call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " +&
      trim(routine_name) + ": must be -1,0,1,2.")
  endif
  if (iflavor==-1) then
    ! FHJ: read flavor directly into iflavor, make sure value is reasonable
    call hdf5_read_int(hidFile, '/mf_header/flavor', iflavor, iError)
    if (iflavor/=1 .and. iflavor/=2) then
      write(sflavor,'(i0)') iflavor
      call die("Illegal flavor = " + TRUNC(sflavor) + " in file " + trim(sFileName))
    endif
  else
    ! FHJ: just check flavor against iflavor2
    call hdf5_require_flavor(hidFile, 'mf_header/flavor', iflavor2, trim(sFileName))
  endif
  call hdf5_require_version(hidFile, 'mf_header/versionnumber', VER_WFN_HDF5, trim(sFileName))
  call h5fclose_f(hidFile, iError)
 
end subroutine read_info
subroutine read_gspace(sFileName,gvec)
  character(len=*), intent(in) :: sFileName
  type(gspace), intent(out) :: gvec
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
  call hdf5_read_int(hidFile, '/mf_header/gspace/ng', gvec%ng, iError)
  call hdf5_read_double(hidFile, '/mf_header/gspace/ecutrho', gvec%ecutrho, iError)
  call hdf5_read_int_array(hidFile, '/mf_header/gspace/FFTgrid', (/3/), gvec%FFTgrid, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine read_gspace
subroutine read_hdf5_gvectors(sFileName, ng, gvec)
  character(len=*), intent(in) :: sFileName
  integer, intent(in) :: ng !< used size of array
  integer, intent(out) :: gvec(:, :) !< (3, ng_bound)
  integer(HID_T) :: hidFile
  integer :: iError
  logical :: bcast_, dont_read_
 
  dont_read_=.false.
  bcast_=.not. dont_read_
  if(peinf%inode == 0) then
    call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
    call hdf5_read_int_array(hidFile, '/mf_header/gspace/components', (/3,ng/), gvec, iError)
    call h5fclose_f(hidFile, iError)
  endif
  if(peinf%npes > 1) then
    if(bcast_) then
      call MPI_BCAST(gvec(1,1), 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
  endif
 
end subroutine read_hdf5_gvectors
subroutine read_symmetry(sFileName,syms)
  character(len=*), intent(in) :: sFileName
  type(symmetry), intent(out) :: syms
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
  call hdf5_read_int(hidFile, '/mf_header/symmetry/ntran', syms%ntran, iError)
  call hdf5_read_int(hidFile, '/mf_header/symmetry/cell_symmetry', syms%cell_symmetry, iError)
  call hdf5_read_int_array(hidFile, '/mf_header/symmetry/mtrx', (/3, 3, 48/), syms%mtrx, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/symmetry/tnp', (/3, 48/), syms%tnp, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine read_symmetry
subroutine read_crystal(sFileName,crys)
  character(len=*), intent(in) :: sFileName
  type(crystal), intent(out) :: crys
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
  call hdf5_read_double(hidFile, '/mf_header/crystal/celvol', crys%celvol, iError)
  call hdf5_read_double(hidFile, '/mf_header/crystal/recvol', crys%recvol, iError)
  call hdf5_read_double(hidFile, '/mf_header/crystal/alat', crys%alat, iError)
  call hdf5_read_double(hidFile, '/mf_header/crystal/blat', crys%blat, iError)
  call hdf5_read_int(hidFile, '/mf_header/crystal/nat', crys%nat, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/crystal/avec', (/3, 3/), crys%avec, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/crystal/bvec', (/3, 3/), crys%bvec, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/crystal/adot', (/3, 3/), crys%adot, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/crystal/bdot', (/3, 3/), crys%bdot, iError)
  allocate(crys%atyp (crys%nat))
  allocate(crys%apos (3,crys%nat))
  call hdf5_read_int_array(hidFile, '/mf_header/crystal/atyp', (/crys%nat/), crys%atyp, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/crystal/apos', (/3, crys%nat/), crys%apos, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine read_crystal
subroutine read_kpoints(sFileName,kp)
  character(len=*), intent(in) :: sFileName
  type(kpoints), intent(out) :: kp
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
  call hdf5_read_int(hidFile, '/mf_header/kpoints/nspin', kp%nspin, iError)
  call hdf5_read_int(hidFile, '/mf_header/kpoints/nspinor', kp%nspinor, iError)
  call hdf5_read_int(hidFile, '/mf_header/kpoints/nrk', kp%nrk, iError)
  call hdf5_read_int(hidFile, '/mf_header/kpoints/mnband', kp%mnband, iError)
  call hdf5_read_int(hidFile, '/mf_header/kpoints/ngkmax', kp%ngkmax, iError)
  call hdf5_read_double(hidFile, '/mf_header/kpoints/ecutwfc', kp%ecutwfc, iError)
  call hdf5_read_int_array(hidFile, '/mf_header/kpoints/kgrid', (/3/), kp%kgrid, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/kpoints/shift', (/3/), kp%shift, iError)
  allocate(kp%ngk (kp%nrk))
  allocate(kp%ifmin (kp%nrk, kp%nspin))
  allocate(kp%ifmax (kp%nrk, kp%nspin))
  allocate(kp%w (kp%nrk))
  allocate(kp%rk (3,kp%nrk))
  allocate(kp%el (kp%mnband,kp%nrk,kp%nspin))
  allocate(kp%occ (kp%mnband,kp%nrk,kp%nspin))
  call hdf5_read_int_array(hidFile, '/mf_header/kpoints/ngk', &
    (/kp%nrk/), kp%ngk, iError)
  call hdf5_read_int_array(hidFile, '/mf_header/kpoints/ifmin', &
    (/kp%nrk, kp%nspin/), kp%ifmin, iError)
  call hdf5_read_int_array(hidFile, '/mf_header/kpoints/ifmax', &
    (/kp%nrk, kp%nspin/), kp%ifmax, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/kpoints/w', &
    (/kp%nrk/), kp%w, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/kpoints/rk', &
    (/3, kp%nrk/), kp%rk, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/kpoints/el', &
    (/kp%mnband, kp%nrk, kp%nspin/), kp%el, iError)
  call hdf5_read_double_array(hidFile, '/mf_header/kpoints/occ', &
    (/kp%mnband, kp%nrk, kp%nspin/), kp%occ, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine read_kpoints
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
! end read/write wfn data
! read/write other
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
subroutine read_hdf5_wfn_gvectors(fname, gvec, ngktot)
  character(len=*) :: fname
  integer, intent(inout) :: gvec(:,:)
  integer, intent(in) :: ngktot
  integer(HID_T) :: file_id
  integer :: error
 
  if(peinf%inode == 0) then
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
    call hdf5_read_int_array(file_id, 'wfns/gvecs', (/3, ngktot/), gvec, error)
    call h5fclose_f(file_id, error)
  endif
  if(peinf%npes > 1) then
    call MPI_Bcast(gvec(1,1), 3*ngktot, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif
 
end subroutine read_hdf5_wfn_gvectors
! end read/write wfn gvectors
! begin read/write wfn data
! end read/write wfn data
! read/write other
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
subroutine read_hdf5_band_complex(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
  character(len=*) :: fname
  complex(DPC), intent(out) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
  integer, intent(in) :: ngk
  integer, intent(in) :: nstot
  integer, intent(in) :: ioffsetk
  integer, intent(in) :: ioffsetb
  real(DP) :: dwfn(2,ngk,nstot,1)
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dataspace_id
  integer(HID_T) :: memspace_id
  integer(HSIZE_T) :: a3(4), offset(4), count(4)
  integer :: error,is
 
  a3(1) = 2
  a3(2) = ngk
  a3(3) = nstot
  a3(4) = 1
  offset(1) = 0
  offset(2) = ioffsetk
  offset(3) = 0
  offset(4) = ioffsetb
  count(1) = 2
  count(2) = ngk
  count(3) = nstot
  count(4) = 1
  if(peinf%inode == 0) then
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
    call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
    CALL h5screate_simple_f(4, count, memspace_id, error)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(dataspace_id, error)
    call h5sclose_f(memspace_id, error)
    call h5fclose_f(file_id, error)
  endif
  wfn(:,:) = cmplx(dwfn(1,:,:,1),dwfn(2,:,:,1),kind=DPC)
 
end subroutine read_hdf5_band_complex
! end read/write wfn data
! read/write other
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
subroutine read_hdf5_band_real(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
  character(len=*) :: fname
  real(DPC), intent(out) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
  integer, intent(in) :: ngk
  integer, intent(in) :: nstot
  integer, intent(in) :: ioffsetk
  integer, intent(in) :: ioffsetb
  real(DP) :: dwfn(1,ngk,nstot,1)
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dataspace_id
  integer(HID_T) :: memspace_id
  integer(HSIZE_T) :: a3(4), offset(4), count(4)
  integer :: error,is
 
  a3(1) = 1
  a3(2) = ngk
  a3(3) = nstot
  a3(4) = 1
  offset(1) = 0
  offset(2) = ioffsetk
  offset(3) = 0
  offset(4) = ioffsetb
  count(1) = 1
  count(2) = ngk
  count(3) = nstot
  count(4) = 1
  if(peinf%inode == 0) then
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
    call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
    CALL h5screate_simple_f(4, count, memspace_id, error)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(dataspace_id, error)
    call h5sclose_f(memspace_id, error)
    call h5fclose_f(file_id, error)
  endif
  wfn(:,:) = dwfn(1,:,:,1)
 
end subroutine read_hdf5_band_real
! end read/write wfn data
! read/write other
!===============================================================================
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
subroutine write_hdf5_header_type(sFileName, sheader, iflavor, kp, gvec, syms, crys)
  character(len=*), intent(in) :: sFileName
  character(len=3), intent(inout) :: sheader
  !> define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(in) :: iflavor
  type(kpoints), intent(in) :: kp
  type(gspace), intent(in) :: gvec
  type(symmetry), intent(in) :: syms
  type(crystal), intent(in) :: crys
  ! set values based on epsilon calculation
  logical :: is_get=.false.
  logical :: wfnflag=.true.
 
  if (peinf%inode == 0) then
    call write_info(TRUNC(sFileName),iflavor)
    call write_kpoints(TRUNC(sFileName),kp)
    call write_gspace(TRUNC(sFileName),gvec)
    call write_symmetry(TRUNC(sFileName),syms)
    call write_crystal(TRUNC(sFileName),crys)
  endif
 
end subroutine write_hdf5_header_type
subroutine write_info(sFileName, iflavor)
  character(len=*), intent(in) :: sFileName
  !> define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(in) :: iflavor
  integer(HID_T) :: hidFile
  integer :: iError, iflavor2
  character(len=16) :: sflavor
  character(len=128) :: routine_name
 
  call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
  iflavor2 = iflavor
  if (iflavor==0) iflavor2 = 1
  routine_name = "write_info"
  ! FHJ: Writing flavor: just write adjusted value in iflavor2
  if (iflavor<0 .or. iflavor>2) then
    write(sflavor,'(i0)') iflavor
    call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " +&
      trim(routine_name) + ": must be 0,1,2.")
  endif
  call hdf5_write_int(hidFile, '/mf_header/versionnumber', VER_WFN_HDF5, iError)
  call hdf5_write_int(hidFile, '/mf_header/flavor', iflavor2, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine write_info
subroutine write_gspace(sFileName,gvec)
  character(len=*), intent(in) :: sFileName
  type(gspace), intent(in) :: gvec
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
  call hdf5_write_int(hidFile, '/mf_header/gspace/ng', gvec%ng, iError)
  call hdf5_write_double(hidFile, '/mf_header/gspace/ecutrho', gvec%ecutrho, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/gspace/FFTgrid', (/3/), gvec%FFTgrid, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine write_gspace
subroutine write_hdf5_gvectors(sFileName, ng, gvec)
  character(len=*), intent(in) :: sFileName
  integer, intent(in) :: ng !< used size of array
  integer, intent(in) :: gvec(:, :) !< (3, ng_bound)
  integer(HID_T) :: hidFile
  integer :: iError
  logical :: bcast_, dont_read_
 
  dont_read_=.false.
  bcast_=.not. dont_read_
  if(peinf%inode == 0) then
    call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/gspace/components', (/3,ng/), gvec, iError)
    call h5fclose_f(hidFile, iError)
  endif
 
end subroutine write_hdf5_gvectors
subroutine write_symmetry(sFileName,syms)
  character(len=*), intent(in) :: sFileName
  type(symmetry), intent(in) :: syms
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
  call hdf5_write_int(hidFile, '/mf_header/symmetry/ntran', syms%ntran, iError)
  call hdf5_write_int(hidFile, '/mf_header/symmetry/cell_symmetry', syms%cell_symmetry, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/symmetry/mtrx', (/3, 3, 48/), syms%mtrx, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/symmetry/tnp', (/3, 48/), syms%tnp, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine write_symmetry
subroutine write_crystal(sFileName,crys)
  character(len=*), intent(in) :: sFileName
  type(crystal), intent(in) :: crys
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
  call hdf5_write_double(hidFile, '/mf_header/crystal/celvol', crys%celvol, iError)
  call hdf5_write_double(hidFile, '/mf_header/crystal/recvol', crys%recvol, iError)
  call hdf5_write_double(hidFile, '/mf_header/crystal/alat', crys%alat, iError)
  call hdf5_write_double(hidFile, '/mf_header/crystal/blat', crys%blat, iError)
  call hdf5_write_int(hidFile, '/mf_header/crystal/nat', crys%nat, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/crystal/avec', (/3, 3/), crys%avec, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/crystal/bvec', (/3, 3/), crys%bvec, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/crystal/adot', (/3, 3/), crys%adot, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/crystal/bdot', (/3, 3/), crys%bdot, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/crystal/atyp', (/crys%nat/), crys%atyp, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/crystal/apos', (/3, crys%nat/), crys%apos, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine write_crystal
subroutine write_kpoints(sFileName,kp)
  character(len=*), intent(in) :: sFileName
  type(kpoints), intent(in) :: kp
  integer(HID_T) :: hidFile
  integer :: iError
 
  call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
  call hdf5_write_int(hidFile, '/mf_header/kpoints/nspin', kp%nspin, iError)
  call hdf5_write_int(hidFile, '/mf_header/kpoints/nspinor', kp%nspinor, iError)
  call hdf5_write_int(hidFile, '/mf_header/kpoints/nrk', kp%nrk, iError)
  call hdf5_write_int(hidFile, '/mf_header/kpoints/mnband', kp%mnband, iError)
  call hdf5_write_int(hidFile, '/mf_header/kpoints/ngkmax', kp%ngkmax, iError)
  call hdf5_write_double(hidFile, '/mf_header/kpoints/ecutwfc', kp%ecutwfc, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/kpoints/kgrid', (/3/), kp%kgrid, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/kpoints/shift', (/3/), kp%shift, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/kpoints/ngk', &
    (/kp%nrk/), kp%ngk, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/kpoints/ifmin', &
    (/kp%nrk, kp%nspin/), kp%ifmin, iError)
  call hdf5_write_int_array(hidFile, '/mf_header/kpoints/ifmax', &
    (/kp%nrk, kp%nspin/), kp%ifmax, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/kpoints/w', &
    (/kp%nrk/), kp%w, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/kpoints/rk', &
    (/3, kp%nrk/), kp%rk, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/kpoints/el', &
    (/kp%mnband, kp%nrk, kp%nspin/), kp%el, iError)
  call hdf5_write_double_array(hidFile, '/mf_header/kpoints/occ', &
    (/kp%mnband, kp%nrk, kp%nspin/), kp%occ, iError)
  call h5fclose_f(hidFile, iError)
 
end subroutine write_kpoints
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
! end read/write wfn data
! read/write other
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
subroutine write_hdf5_wfn_gvectors(fname, gvec, ngktot)
  character(len=*) :: fname
  integer, intent(inout) :: gvec(:,:)
  integer, intent(in) :: ngktot
  integer(HID_T) :: file_id
  integer :: error
 
  if(peinf%inode == 0) then
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    call hdf5_write_int_array(file_id, 'wfns/gvecs', (/3, ngktot/), gvec, error)
    call h5fclose_f(file_id, error)
  endif
 
end subroutine write_hdf5_wfn_gvectors
! end read/write wfn gvectors
! begin read/write wfn data
! end read/write wfn data
! read/write other
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
subroutine write_hdf5_band_complex(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
  character(len=*) :: fname
  complex(DPC), intent(in) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
  integer, intent(in) :: ngk
  integer, intent(in) :: nstot
  integer, intent(in) :: ioffsetk
  integer, intent(in) :: ioffsetb
  real(DP) :: dwfn(2,ngk,nstot,1)
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dataspace_id
  integer(HID_T) :: memspace_id
  integer(HSIZE_T) :: a3(4), offset(4), count(4)
  integer :: error,is
 
  a3(1) = 2
  a3(2) = ngk
  a3(3) = nstot
  a3(4) = 1
  offset(1) = 0
  offset(2) = ioffsetk
  offset(3) = 0
  offset(4) = ioffsetb
  count(1) = 2
  count(2) = ngk
  count(3) = nstot
  count(4) = 1
  dwfn(1,:,:,1) = real(wfn(:,:))
  dwfn(2,:,:,1) = aimag(wfn(:,:))
  if(peinf%inode == 0) then
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
    CALL h5screate_simple_f(4, count, memspace_id, error)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(dataspace_id, error)
    call h5sclose_f(memspace_id, error)
    call h5fclose_f(file_id, error)
  endif
 
end subroutine write_hdf5_band_complex
! end read/write wfn data
! read/write other
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
subroutine write_hdf5_band_real(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
  character(len=*) :: fname
  real(DPC), intent(in) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
  integer, intent(in) :: ngk
  integer, intent(in) :: nstot
  integer, intent(in) :: ioffsetk
  integer, intent(in) :: ioffsetb
  real(DP) :: dwfn(1,ngk,nstot,1)
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dataspace_id
  integer(HID_T) :: memspace_id
  integer(HSIZE_T) :: a3(4), offset(4), count(4)
  integer :: error,is
 
  a3(1) = 1
  a3(2) = ngk
  a3(3) = nstot
  a3(4) = 1
  offset(1) = 0
  offset(2) = ioffsetk
  offset(3) = 0
  offset(4) = ioffsetb
  count(1) = 1
  count(2) = ngk
  count(3) = nstot
  count(4) = 1
  dwfn(1,:,:,1) = wfn(:,:)
  if(peinf%inode == 0) then
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
    CALL h5screate_simple_f(4, count, memspace_id, error)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
    call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
    call h5dclose_f(dataset_id, error)
    call h5sclose_f(dataspace_id, error)
    call h5sclose_f(memspace_id, error)
    call h5fclose_f(file_id, error)
  endif
 
end subroutine write_hdf5_band_real
! end read/write wfn data
! read/write other
!
!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================
! begin read/write header
! end read/write header
! begin read/write wfn gvectors
! end read/write wfn gvectors
! begin read/write wfn data
! end read/write wfn data
subroutine read_hdf5_bands_block(file_id, kp, nbownmax, nbownactual, does_it_ownb, ib_first, wfnsout, ioffset)
  integer(HID_T), intent(in) :: file_id
  type(kpoints), intent(in) :: kp
  integer, intent(in) :: nbownmax
  integer, intent(in) :: nbownactual !< how many bands I own
  logical, intent(in) :: does_it_ownb(:,:)
  integer, intent(in) :: ib_first !< first band of the set of bands I own
  real(DP), intent(out) :: wfnsout(:,:,:) !< (SUM(kp%ngk), kp%nspin*kp%nspinor, nbownactual)
  integer, optional, intent(in) :: ioffset
  real(DP), allocatable :: wfndata(:,:,:,:)
  integer(HID_T) :: plist_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: dataspace
  integer(HID_T) :: memspace
  integer(HSIZE_T) :: count(4), offset(4)
  integer :: error
  integer :: ipe, reader
  integer :: nread
  integer :: ngktot
  integer :: ioffset_
  integer, allocatable :: ranks(:)
  integer :: icount
  integer :: mpiworldgrp, mpigrp, bandcomm, ib, ib_, max_bands_read, bands_read, bands_read_max, is
  integer, parameter :: max_bytes_read = 536870912 ! don`t read/send more than 1/2 GB at a time
  ! 0=native 1; 1=manual group comm; 2=manual send/recvs
  integer, parameter :: comm_style=0
  logical :: do_read
  !real(DP) :: mem
  !integer :: nmpinode
 
  call logit('Reading HDF5 WFNs by blocks')
  allocate(ranks (peinf%npes))
  ioffset_ = 0
  if(present(ioffset)) then
    ioffset_ = ioffset
  endif
  ngktot = SUM(kp%ngk)
  nread = peinf%npools
  call h5dopen_f(file_id, 'wfns/coeffs', dset_id, error)
  call h5dget_space_f(dset_id, dataspace, error)
  ! find lowest rank PE that owns the first band of this group of bands
  reader = -1
  if (nbownactual>0) then
    do ipe = 1, peinf%npes
      if(does_it_ownb(ib_first,ipe)) then
        reader = ipe - 1
        exit
      endif
    enddo
    if (reader==-1) call die("Cannot find reader in read_hdf5_bands_block", only_root_writes=.true.)
    if (comm_style==1) then
      ! We use MPI_Bcast with 1 groups
      icount = 0
      do ipe = 1, peinf%npes
        if(does_it_ownb(ib_first,ipe)) then
            icount = icount + 1
            ranks(icount) = ipe - 1
        endif
      enddo
      call MPI_Comm_Group(MPI_COMM_WORLD, mpiworldgrp, mpierr)
      call MPI_Group_Incl(mpiworldgrp, icount, ranks, mpigrp, mpierr)
      call MPI_Comm_Create(MPI_COMM_WORLD, mpigrp, bandcomm, mpierr)
    endif
  else
    if (comm_style==1) then
      ! FHJ: Note that MPI_Comm_Create must be called by everyone in MPI_COMM_WORLD!
      call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, bandcomm, mpierr)
    endif
  endif
  ! FHJ: read at most max_bands_read bands to avoid 1/1 buffer overflow.
  ! Here, 1*kp%nspin*kp%nspinor*dble(ngktot)*8d0 is the size of a
  ! single band, including G-vectors from all k-points.
  max_bands_read = min(nbownmax, &
    int(max_bytes_read/(1*kp%nspin*kp%nspinor*dble(ngktot)*8d0)))
  if (max_bands_read==0) then
    max_bands_read = 1
    if (peinf%inode==0) then
      write(0,*)
      write(0,'(a)') 'WARNING: could not honor limit of '
      write(0,'(f0.3,a)') max_bytes_read/1024d0/1024d0,' MB per chunk when'
      write(0,'(a)') 'reading HDF5 WFN file. Using chunks of '
      write(0,'(f0.3,a)') (kp%nspin*kp%nspinor*1*dble(ngktot)*8d0)/1024d0/1024d0,' MB.'
      write(0,*)
    endif
  endif
  !write(6,*) 'max_bands_read', max_bands_read
  allocate(wfndata (1, ngktot, kp%nspin*kp%nspinor, max_bands_read))
  ib = 1

  ! OAH this is likely to be where we make modifications to the read-ins
  do while (ib<=nbownmax)
    bands_read = max(min(nbownactual, ib-1 + max_bands_read) - ib + 1, 0)
    bands_read_max = max(min(nbownmax, ib-1 + max_bands_read) - ib + 1, 0)
    count(1) = 1
    count(2) = ngktot
    count(3) = kp%nspin*kp%nspinor
    count(4) = bands_read
    call h5screate_simple_f(4, count, memspace, error)
    do_read = bands_read>0.and.(peinf%inode==reader.or.comm_style==0)
    if (do_read) then
      offset(1) = 0
      offset(2) = 0
      offset(3) = 0
      offset(4) = (ib_first-1)+ioffset_+(ib-1)
      if (peinf%verb_debug .and. peinf%inode==reader) then
        write(6,'(4(a,i0,1x))') 'ib=',ib,'bands_read=',bands_read,'offset=',offset(3),'ib_first=',ib_first
      endif
    else
      offset(:) = 0
      call H5sselect_none_f(memspace,error)
    endif
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (.not.do_read) then
      call H5sselect_none_f(dataspace,error)
    endif
    if (peinf%verb_debug .and. peinf%inode==reader) then
      write(6,'(a,i0,a)') 'ib=',ib,' before read!'
      !call procmem(mem,nmpinode)
      !write(6,'(a,f0.3,a)') 'Memory available: ', mem/(1024d0**2),' MB per MPI rank.'
    endif
    !if (peinf%inode==reader) write(6,'(a)') '[1]'
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    !if (peinf%inode==reader) write(6,'(a)') '[2]'
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !if (peinf%inode==reader) write(6,'(a)') '[3]'
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata(:,:,:,:), count, error, memspace, dataspace, xfer_prp=plist_id)
    !if (peinf%inode==reader) write(6,'(a)') '[4]'
    call h5pclose_f(plist_id, error)
    !if (peinf%inode==reader) write(6,'(a)') '[5]'
    if (peinf%verb_debug .and. peinf%inode==reader) then
      write(6,'(a,i0,a)') 'ib=',ib,' after read!'
    endif
    call h5sclose_f(memspace, error)
    if (do_read) then
      do is = 1, kp%nspin*kp%nspinor
        wfnsout(:,is,ib:ib+bands_read-1) = &
          real(wfndata(1,:,is,1:bands_read),kind=DP)
      enddo
    endif
    if (bands_read>0) then
      ! FHJ: No manual distribution is necessary for comm_style==0
      if (comm_style>0) call logit('Sending bands')
      if (comm_style==2) then
        if (peinf%inode==reader) then
          do ipe = 1, peinf%npes
            if(does_it_ownb(ib_first,ipe) .and. ipe-1 .ne. peinf%inode) then
              call MPI_Send(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, &
                MPI_DOUBLE_PRECISION, ipe-1, 0, MPI_COMM_WORLD, mpierr)
            endif
          enddo
        else
          call MPI_Recv(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, &
            MPI_DOUBLE_PRECISION, reader, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
        endif
      elseif (comm_style==1) then
        call MPI_Bcast(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, &
          MPI_DOUBLE_PRECISION, 0, bandcomm, mpierr)
      endif
    endif
    ib = ib + bands_read_max ! so this is like a one-shot thing for each task
  enddo
  if (comm_style==1.and.nbownactual>0) then
    call MPI_Comm_free(bandcomm, mpierr)
    call MPI_Group_free(mpigrp, mpierr)
  endif
  if(allocated(ranks))then;deallocate(ranks);endif
  if(allocated(wfndata))then;deallocate(wfndata);endif
  call h5sclose_f(dataspace, error)
  call h5dclose_f(dset_id, error)
 
end subroutine read_hdf5_bands_block
! read/write other
!===============================================================================
  !> Create the appropriate structures that hold info about the mean-field
  !! calculation in an 1 file. No data is actually written.
  subroutine setup_hdf5_mf_file(fname, create_file)
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: create_file
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer :: error
    logical :: create_file_
   
    create_file_ = .true.
    if (present(create_file)) create_file_ = create_file
    if (create_file_) then
      call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
    else
      call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    endif
    call h5gcreate_f(file_id, '/mf_header', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/gspace', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/symmetry', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/crystal', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/kpoints', group_id, error)
    call h5gclose_f(group_id, error)
    call h5fclose_f(file_id, error)
   
  end subroutine setup_hdf5_mf_file
  !> Create the appropriate structures that hold info about the WFN
  !! coefficients in an 1 file. No data is actually written.
  subroutine setup_hdf5_wfn_file(fname, iflavor, kp)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: dataset_id
    integer(HSIZE_T) :: a3(4)
    integer :: error
   
    call setup_hdf5_mf_file(fname)
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    call h5gcreate_f(file_id, '/wfns', group_id, error)
    call h5gclose_f(group_id, error)
    if (iflavor/=1 .and. iflavor/=2) then
      write(0,*) 'ERROR: got iflavor=', iflavor
      call die('Internal error: invalid flavor in setup_hdf5_wfn_file.', only_root_writes=.true.)
    endif
    a3(1) = iflavor
    a3(2) = sum(kp%ngk)
    a3(3) = kp%nspin*kp%nspinor
    a3(4) = kp%mnband
    CALL h5screate_simple_f(4, a3, dataspace_id, error)
    call h5dcreate_f(file_id, '/wfns/coeffs', H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
    CALL h5sclose_f(dataspace_id, error)
    call h5fclose_f(file_id, error)
   
  end subroutine setup_hdf5_wfn_file
  !> A high-level wrapper for write_*_header* functions
  subroutine write_hdf5_mf_header(fname, mf)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(in) :: mf
    character(len=3) :: sheader
   
    sheader = mf%sheader
    call write_hdf5_header_type(fname, sheader, mf%iflavor, &
      !mf%kp, mf%gvec, mf%syms, mf%crys, version=mf%version)
      mf%kp, mf%gvec, mf%syms, mf%crys)
   
  end subroutine write_hdf5_mf_header
  !> A high-level wrapper for read_*_header* functions
  !! Note: the optional fields are ignored for now.
  subroutine read_hdf5_mf_header(fname, mf, iflavor, sheader, warn, dont_warn_kgrid)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(out) :: mf
    integer, intent(in), optional :: iflavor
    character(len=3), intent(in), optional :: sheader
    logical, intent(in), optional :: warn
    logical, intent(in), optional :: dont_warn_kgrid
   
    if (present(sheader)) then
      mf%sheader = sheader
    else
      mf%sheader = 'GET'
    endif
    if (present(iflavor)) then
      mf%iflavor = iflavor
    else
      mf%iflavor = -1
    endif
    call read_hdf5_header_type(fname, mf%sheader, mf%iflavor, &
      !mf%kp, mf%gvec, mf%syms, mf%crys, version=mf%version, sdate=mf%sdate, stime=mf%stime)
      mf%kp, mf%gvec, mf%syms, mf%crys)
    !FHJ: FIXME - Implement in WFN file
    mf%sheader = 'WFN'
    mf%stime = 'N/A'
    mf%sdate = 'N/A'
   
  end subroutine read_hdf5_mf_header
end module wfn_io_hdf5_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
