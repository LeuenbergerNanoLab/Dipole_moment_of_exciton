

 Module Module_XSF
        implicit none
        real(8), allocatable        :: Pot(:)                    ! excitonic wave function
        real(8), allocatable        :: Pos(:,:)                  ! grid mesh
        integer                     :: Nxpot,Nypot,Nzpot         ! mesh size
        integer                     :: Npot                      ! total number of mesh points
        character(30)               :: name 
        real(8)                     :: Point(3)
        real(8)                     :: Period(3,3)               ! lattice vectors
        integer                     :: Nat                       ! number of atoms
        character(2),allocatable    :: Atoms(:)                  ! atom names
        real(8)     ,allocatable    :: Coord(:,:)                ! positions of the atoms
        real(8)                     :: dx(3),dy(3),dz(3)
        logical                     :: lPot
        real(8)                     :: hole(3)                   ! position of the hole
        real(8)                     :: dipole(3)                 ! dipole moment of exciton
 Contains



  subroutine read_xsf
    integer            :: i,j,k
    integer            :: jj
    print *,'read_xsf: ',name
    lPot = .true.
    open(unit=2,file=trim(adjustl(name)))
    do i=1,46                                   ! xsf created by Yambo has 46 comment lines
     read(2,*)
    enddo
    read(2,*)
    read(2,*) 
    read(2,*) Period(1:3,1)
    read(2,*) Period(1:3,2)
    read(2,*) Period(1:3,3)
    read(2,*) 
    read(2,*) Nat
    allocate(Coord(3,Nat))
    allocate(Atoms(Nat))
    do i=1,Nat
     read(2,*) Atoms(i),Coord(1:3,i)
    enddo
    hole(1:3) = Coord(1:3,1)                      ! first atom = hole
    print *,'hole=',hole
    if(lPot) then
     read(2,*)
     read(2,*)
     read(2,*)
     read(2,*) Nxpot,Nypot,Nzpot
     Npot = Nxpot*Nypot*Nzpot
     read(2,*) Point
     Point = -hole                      ! make hole = (0,0,0)
     read(2,*) Period(1:3,1)
     read(2,*) Period(1:3,2)
     read(2,*) Period(1:3,3)
     dx(1:3) = Period(1:3,1)/dfloat(Nxpot-1)
     dy(1:3) = Period(1:3,2)/dfloat(Nypot-1)
     dz(1:3) = Period(1:3,3)/dfloat(Nzpot-1)
     allocate(Pot(Npot))
     allocate(Pos(3,Npot))
     print *,'Pot: allocated ',allocated(Pot)
     print *,'Pos: allocated ',allocated(Pos)
     read(2,*) (((Pot(i+(j-1)*Nxpot+(k-1)*Nxpot*Nypot),i=1,Nxpot),j=1,Nypot),k=1,Nzpot)
     print *,'Nxpot=',Nxpot
     print *,'Nypot=',Nypot
     print *,'Nzpot=',Nzpot
     print *,'Npot=',Npot
     dx(1:3) = Period(1:3,1)/dfloat(Nxpot-1)
     dy(1:3) = Period(1:3,2)/dfloat(Nypot-1)
     dz(1:3) = Period(1:3,3)/dfloat(Nzpot-1)
     print *,'dx(1:3)=',dx(1:3)
     print *,'dy(1:3)=',dy(1:3)
     print *,'dz(1:3)=',dz(1:3)
     jj = 0
     do k=1,Nzpot                                          ! calculate the real space mesh
      do j=1,Nypot
       do i=1,Nxpot
        jj = jj + 1
        if(jj > Npot) then
         print *,'read_potential: error jj > Npot'
         stop
        endif
        Pos(1:3,jj) = (i-1)*dx(1:3) + (j-1)*dy(1:3) + (k-1)*dz(1:3) + Point(1:3)
       enddo
      enddo
     enddo
    endif        ! lPot
    close(unit=2)
    print *,'read_xsf:  finish'
6   format(6E14.6)
11  format(I3,2x,A2,3F17.5) 
  end subroutine read_xsf



  subroutine calc_dipole_moment
   integer       :: j
   real(8)       :: r(3)                   ! position operator
   real(8)       :: nWF                    ! norm
   dipole = 0.d0
   nWF = 0.d0
   do j=1,Npot
    nWF = nWF + Pot(j)
   enddo
   nWF = nWF/dfloat(Npot)
   print *,'norm=',nWF
   do j=1,Npot
    r = (/ pos(1,j), pos(2,j), pos(3,j) /)
    dipole = dipole + r*Pot(j)
   enddo
   dipole = dipole/dfloat(Npot)/nWF
   print 1,dipole
1  format('dipole(1)=',E12.5/ &
         ' dipole(2)=',E12.5/ &
         ' dipole(3)=',E12.5)   
  end subroutine calc_dipole_moment



 end module Module_XSF



