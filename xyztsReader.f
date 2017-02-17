c!----------------------------------------------------------------------
        module xyztsInp

        integer idomain                         !The domain type
        integer nbl
        real*8  yplus, bl_thick, ratio
        real*8  pitch, rpin, rmax, rcentre
        real*8,allocatable:: probe(:,:)
        integer N, S
        real*8  axial
        real*8  yc, zc          !center of cross plane

        integer nhom                            !for pipe

        end module
c!----------------------------------------------------------------------
        subroutine xyztsReader()
c!----------------------------------------------------------------------
c!
c!...This subroutine is used to read the input
c!
c! Jun Fang                                             Fall, 2014
c!----------------------------------------------------------------------
        use xyztsInp
        integer ierror
        integer i,j

        open(unit = 10, file = 'xyztsGen.inp', iostat = ierror)
        if (ierror .ne. 0) stop"Error opening input file"

        read(10,*)
        read(10,*) idomain
        read(10,*)
        if (idomain .eq. 1) then        !read subchannel parameters
           write(*,*) 'The domain type: Subchannel'
           read(10,*) N
           write(*,*) 'Numbers of centres = ', N
           read(10,*) S
           write(*,*) 'Number of subchannels = ', S
           read(10,*)
           read(10,*) yplus, bl_thick, nbl
           write(*,*) 'y+ = ', yplus
           write(*,*) 'Thickness of BL region = ', bl_thick
           write(*,*) 'Number of boundary layers = ', nbl
           read(10,*) ratio
           write(*,*) 'BL growth ratio = ', ratio
           read(10,*) pitch
           write(*,*) 'Pitch = ', pitch
           read(10,*) rpin
           write(*,*) 'Fuel rod radius = ', rpin
           read(10,*) axial
           write(*,'(A, E12.3)') ' X position of probe plane: ', axial
           read(10,*)
           read(10,*) yc, zc
 
           rcentre      = ((pitch/2.0)**2.0+(pitch/2.0)**2.0)**0.5


        endif

        if (idomain .eq. 2) then        !read pipe parameters
           write(*,*) 'The domain type: Pipe'
           do i = 4,17                  !skip suchannel parameters
              read(10,*)
           enddo
           read(10,*) rmax
           write(*,*) 'Pipe radius = ', rmax
           read(10,*) nhom
           write(*,*) 'number of homogeneous points =', nhom
           read(10,*) axial
           write(*,'(A, E12.3)') 'X position of probe plane: ', axial
        endif

        end !of the subroutine

c!----------------------------------------------------------------------
