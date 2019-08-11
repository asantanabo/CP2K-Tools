program cubeaverage


integer :: i,j,k,m, nx,ny,nz, narg, ln, err, natms
real(8), dimension(:,:,:), allocatable :: phi3d, phi3d_0
real(8), dimension(:), allocatable :: x,y,z
real(8), dimension(3) :: or, dl
real(8), parameter :: au=0.529177
real(8) :: tmp
character(256) :: filebox,filex,filey,filez,filename,refname
character(80) :: buffer
logical :: xav,yav,zav

narg=command_argument_count()

if (.not.(narg.eq.2)) then
 write(*,*) 'usage:'
 write(*,*) 'cubeaverage -xy pot_file '
 stop 
endif
 

call get_command(filebox, ln, err)
if (ln<0 .or. err>0) then
  stop 'Internal error: command line too long'
endif

xav = .true.
yav = .true.
zav = .true.

k = index(filebox,"-xy")
if (k > 0) then
  zav = .false.  
endif
k = index(filebox,"-xz")
if (k > 0) then
  yav = .false.  
endif
k = index(filebox,"-yz")
if (k > 0) then
  xav = .false.  
endif

call get_command_argument(2,filename,ln,err)

open(110,file=trim(filename))
read(110,*) 
read(110,*) 
read(110,*) natms
read(110,*) nx, dl(1) 
read(110,*) ny, tmp, dl(2) 
read(110,*) nz, tmp, tmp, dl(3) 

!print*, natms, nx, ny, nz
!print*, dl
allocate(phi3d(nx,ny,nz))


do i = 1, natms
  read(110,*)
end do


do i = 1, nx
  do j = 1, ny
     read(110,'(6E13.5)') (phi3d(i,j,k), k=1,nz)
  end do
end do

close(110)

!integrate part
if (.not.xav) then 
  allocate(x(nx))
  x = 0.d0
  do i=1,nx
    do j=1,ny
      do k=1,nz
         x(i) = x(i) + phi3d(i,j,k) 
      enddo
    enddo
    write(*,*) i*dl(1), x(i)/(ny*nz)
  enddo
endif

if (.not.yav) then 
  allocate(x(ny))
  x = 0.0d0
  do j=1,ny
    do i=1,nx
      do k=1,nz
         x(j) = x(j) + phi3d(i,j,k) 
      enddo
    enddo
    write(*,*) j*dl(2), x(j)/(nx*nz)
  enddo
endif

if (.not.zav) then 
  allocate(x(nz))
  x = 0.0d0
  do k=1,nz
    do i=1,nx
      do j=1,ny
         x(k) = x(k) + phi3d(i,j,k) 
      enddo
    enddo
    write(*,*) k*dl(3), x(k)/(nx*ny)
  enddo
endif

deallocate(x)
deallocate(phi3d)

end program cubeaverage
