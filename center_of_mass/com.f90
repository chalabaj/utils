! this program calculate center of mass for each geometrz in input file
! there are 3 arguments given bz bash script (input file, number of atoms in geometry and number of geometries in inputfile]
! calculated center of mas is shifted to the origin 0,0,0 -> all coordinates are shifted bz center of mass position.

program center_of_mass
  implicit none
	INTEGER Natoms,Ngeoms,Nlines,i_at,i_geom
  
  REAL*8                        :: Mtot,Xcom,Ycom,Zcom,mx,my,mz,dist,distcom
	REAL*8, allocatable           :: x(:),y(:),z(:),Mass(:)  
  CHARACTER(len=2),allocatable  :: Name(:)
  CHARACTER(20)                 :: arg1,arg2,arg3,inputfile,outputfile
!------------------------------------------------------------------------------
  
  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)  
	read(arg1,'(A)')inputfile
 	read(arg2,'(I10)')Natoms
 	read(arg3,'(I10)')Nlines
  Ngeoms=Nlines/(Natoms+2)
  print *,'Ngeoms, Nlines: ',Ngeoms,Nlines
    
  inputfile=trim(inputfile) 
  print *,'Input file: ',inputfile  
  outputfile='com-'//inputfile
  print *,'COM geometries will be saved to: ',outputfile
  allocate ( x(natoms)) 
  allocate ( y(natoms)) 
  allocate ( z(natoms)) 
  allocate ( Name(natoms)) 
  allocate ( Mass(natoms))   

  open(101,file=inputfile,status='OLD') !open file .xyz
  open(102,file=outputfile,status='REPLACE')  
  print *,
  do i_geom=1,Ngeoms
    Mtot=0
    Xcom=0
    mx=0
    Ycom=0
    my=0
    Zcom=0
    mz=0
    read(101,*)
    read(101,*)
    do i_at=1,natoms 
      read(101,*)Name(i_at),x(i_at),y(i_at),z(i_at)
      Mass(i_at)=masses(Name(i_at))
      Mtot=Mtot+Mass(i_at)

      mx=mx + (x(i_at)*Mass(i_at))  !sum of x*m
      my=my + (y(i_at)*Mass(i_at))
      mz=mz + (z(i_at)*Mass(i_at))

    end do
    dist=sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2 + (z(1)-z(2))**2) 
    Xcom=mx/Mtot
    Ycom=my/Mtot
    Zcom=mz/Mtot
!    write(*,30)'mtot=',Mtot, mx,my,mz,Xcom,Ycom, Zcom
   
   

    write(102,*)natoms
    write(102,1)'Center of mass position: [',Xcom,' ',Ycom,' ',Zcom,']'
1  format(A5,ES16.5,A2,ES16.5,A2,ES16.5,A2)    
! move COM to 0,0,0 -> shift each atomwith respect to shift of the center of mass 
    do i_at=1,natoms   
      x(i_at)=x(i_at)-Xcom
      y(i_at)=y(i_at)-Ycom
      z(i_at)=z(i_at)-Zcom
      write(102,2)Name(i_at),x(i_at),y(i_at),z(i_at)
2     format(1A2,3E16.8)   
    end do  
  distcom=sqrt((x(1)-x(2))**2 + (y(1)-y(2))**2 + (z(1)-z(2))**2) 
  write(*,3)dist,distcom
3 format(2E16.5)
  end do         
 close(101)
 close(102)
        
CONTAINS         

  REAL*8 FUNCTION masses(Name)
  IMPLICIT NONE
  Character(len=2),INTENT(IN) :: Name
  REAL*8 :: mass
  SELECT case (trim(Name))
   case ('O')
    mass=15.9994D0
   case ('H')
    mass=1.0079D0
  end SELECT
  masses=mass 
 END FUNCTION masses
 
 
end program center_of_mass      