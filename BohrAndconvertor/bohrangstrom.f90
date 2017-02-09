 program bta_convertor
 IMPLICIT NONE
 REAL*8           :: x,y,z
 INTEGER          :: Ngeoms,j,Natoms,i_atom,Nargs
 CHARACTER(50)    :: inputfile(1000),outputfile(1000),arg
 CHARACTER(100)   :: description
 CHARACTER(2)     :: name,units
 REAL*8,parameter :: bTa=0.52917720859

 ! Error handling-----------------------------------------------------------------
 Nargs=command_argument_count()
 if ( Nargs.LT.2 )then
  call Print_help(1)
 end if 
 
 units = '-N'   ! which unit to convert    
 call get_command_argument(1, units)    
 if ( units.EQ.'-N' ) then
  call Print_help(2)
 end if 

 !Reading input files-------------------------------------------------------------
 Ngeoms = Nargs - 1  
 j=1               ! j-th geometry
      do while (j .LE. Ngeoms)
          call get_command_argument(j+1, arg)     ! j+1 because first argument is controling the units and the rest are the files to process  
	        read(arg,'(A)')inputfile(j)
          
          if ( units.EQ.'-B' )then
              outputfile(j)=adjustl('Ang-'//inputfile(j))
          else if ( units.EQ.'-A' ) then
              outputfile(j)=adjustl('Bohr-'//inputfile(j))
          end if
          j=j+1
      end do

!Geometry-------------------------------------------------------------
  do j=1,Ngeoms,1  
      open(110,file=inputfile(j),status='OLD') 
      read(110,*)Natoms
      read(110,*)description
      
      open(111,file=outputfile(j),status='REPLACE')   
      write(111,*)Natoms
      write(111,*)trim(description)   
          
      do i_atom=1,Natoms,1
        read(110,*)name,x,y,z

        if ( units.EQ.'-B' )then
        x = x * bTa
        z = z * bTa
        y = y * bTa
        else if ( units.EQ.'-A' ) then
        x = x / bTa
        z = z / bTa
        y = y / bTa        
        end if 
               
        write(111,2)name,x,y,z ! writing XYZ to movie file   
2       format(1A,3E16.8)   
     
      end do
      close(110)
      close(111)
            
  end do
  
 end program bta_convertor
 
 SUBROUTINE Print_help(err)
 IMPLICIT  NONE
 integer :: err
 
 if (err.EQ.1) then
          write(*,*)'Not enough parameters.'
          write(*,*)'Usage ./bta -B file file file..... to convert from Bohr to Angstroem units'
          write(*,*)'Usage ./bta -A file file file..... to convert from Angstroem to Bohr units'
 else if (err.EQ.2) then
          write(*,*)'Units of original file are not set.'
          write(*,*)'Usage ./bta -B file file file..... to convert from Bohr to Angstroem units'
          write(*,*)'Usage ./bta -A file file file..... to convert from Angstroem to Bohr units'
 end if 
 STOP          
 END SUBROUTINE Print_help