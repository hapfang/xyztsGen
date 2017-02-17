        Program XYZTS_generator
!
!----------------------------------------------------------------------
!
!       This routine is developed to generate probes which can be used
!       to collect the data during PHASTA simulations for post-
!       statistical analysis (Jun Fang, Feb 2017)
!
!       Note: x direction is the vertical or streamwise direction for
!       the current subchannel studies
!
!       List of variables
!       npr             :# of probes around a fuel rod (1st layer)
!       nly             :total # of probe layers
!       npb             :total # of probes 
!----------------------------------------------------------------------
        use xyztsInp
        implicit none

        real*8  r, r1, r2, dr
        real*8  theta, theta1, dtheta, dtheta1
        real*8  x, y, z, y0, z0, pi
        real*8  Angles(1:10, 1:256)   ! Angles for different regions
        real*8  Radii(1:256)

        integer i, j, k, npr, npr1, npb, nly, ily
        integer n1, n2 
        integer Nar, Nrad, Nlayer1
        integer Nregion(1:20), i_max(20)
       
        logical istrue

!       Load the input parameters
        call    xyztsReader()

        pi      = 4.0*atan(1.0)
        npr     = 40
        npr1    = npr/4         !# of probes along 1/4 rod perimeter            

        allocate(probe(npr,3))
!       Open data file to save the probe coordinates
        open(13,file="xyzts.dat")

        if (idomain.eq.1) then

        ily     = 0
        r       = rpin + yplus          !position of 1st layer probe
        r1      = rcentre - 140.0*yplus !below this is central region
        dr      = yplus
        y0      = pitch/2.0d0
        z0      = pitch/2.0d0

        do while (r .lt. r1)            !start radial distance loop
           probe  = 0.0d0
           theta1 = pi                  !beginning angle
           theta  = pi/2.0              !available view-angle
           dtheta = pi/2.0/real(npr1+1) !increment angle

           istrue = .false.                        
           do while(.not.istrue) !determine the view angle and beginning angle
              y = y0 + r*cos(theta1) 
              z = z0 + r*sin(theta1)
              if((y.lt.yc).or.(z.lt.zc)) then
                theta1  = theta1 + 0.1*pi/90.0  !update beginning angle
                theta   = theta  - 0.2*pi/90.0  !update/reduce view-angle
                if(theta.lt.0.0d0) exit 
              elseif(y.gt.yc .and. z.gt.zc) then
                istrue  = .true.                !point within limits
              endif
            enddo 

            dtheta1 = theta/real(npr1 + 1)  !update increment angle
            if(theta.lt.0.0d0) exit
            ily         = ily +1
            probe(:,1)  = axial
            do k = 1, npr1 ! start along circle sectors
               probe(k,2) = y0 + r*cos(theta1 + dtheta1)
               probe(k,3) = z0 + r*sin(theta1 + dtheta1)               
               theta1     = theta1 + dtheta1
            enddo ! end loop along sectors
            probe((npr1+1):(2*npr1),2)   = 2.0*yc - probe(1:npr1,2)    
            probe((npr1+1):(2*npr1),3)   = probe(1:npr1,3)
            probe((2*npr1+1):(3*npr1),2) = 2.0*yc - probe(1:npr1,2)
            probe((2*npr1+1):(3*npr1),3) = 2.0*zc - probe(1:npr1,3)
            probe((3*npr1+1):(4*npr1),2) = probe(1:npr1,2)
            probe((3*npr1+1):(4*npr1),3) = 2.0*zc - probe(1:npr1,3)
            call write_probe(npr)

            if(r .le. (rpin+bl_thick)) then
               dr = dr*ratio
            elseif(r .le. (rpin+rcentre)/2.0) then
               dr = 30.0*yplus
            else
               dr = 40.0*yplus
            endif
            r     = r + dr
        enddo ! end while loop

        probe(:,:)  = 0.0d0 
        probe(:,1)  = axial
        n2 = 0
        do i = 1, 4
           n1          = int(npr*(i-1)*3/20)
           if(i==1) n1 = int(npr*i/10)
           if(i==4) n1 = npr - n2 
           dtheta      = 2.0*pi/real(n1)
           r2          = real(i)*35.0*yplus
           do j = 1, n1
              probe(n2+j,2) = yc + r2*cos(real(j-0.5)*dtheta)
              probe(n2+j,3) = zc + r2*sin(real(j-0.5)*dtheta)
           enddo
           n2   = n2 + n1
        enddo

        call write_probe(npr)
        close(13)
        nly     = ily + 1
        npb     = nly * npr1 * N
        Print *, "radial layers", nly
        write(*,*) 'Total number of probes: ', npb
        Print * , "N,s",N,s, rcentre-rpin

        end if !if loop of domain type


        if (idomain.eq.2) then  ! Pipe option here

! Create an array for radial locations
        
! Number of regions:
         Nar = 3
           
! Create an array for theta angles for all regions:
         do j = 1, Nar
          i_max(j) = int(real(nhom)/2.0**real(j-1))
          do i = 1, i_max(j)
            Angles(j,i) = real(i-1)/real(i_max(j))*2.0*pi           
!           write(*,*) j, i, Angles(j,i)
          end do
         end do

! Radial locations: (ignore BL for now!)
! For uniform in R spacing and preserving total homogeneous directions
! we will need to have a geometric progression in terms of number of
! layers in each "angle" region
         
! First let's compute number of the radius coordinates
        Nrad = 0
! Number of layers for the dense near wall region (will be spaced as BL
! in the future):
        Nlayer1 = 4
        do j = 1, Nar
         Nregion(j) = Nlayer1*int(2.0**(j-1))
         Nrad = Nrad + Nregion(j)
        end do
        write(*,*) 'Nrad = ', Nrad

! Split and fill with coordinates:
        do j = 1, Nrad
!         Radii(j) = rmax*(1.0 - (real(j)-0.5)/real(Nrad))
         Radii(j) = rmax*cos((real(j)-0.5)/real(Nrad)*0.5*pi)
        end do

! See if we can print out the final result:        
! (check the sequence later)
        Nrad = 0
        do j = 1, Nar   ! Loop over the regions
         do k = 1, Nregion(j)  ! Loop over the radial positions within the region
          do i = 1, i_max(j)   ! Loop over the tangential angles
           write(13,*) axial, Radii(Nrad+k)*cos(Angles(j,i))
     &                      , Radii(Nrad+k)*sin(Angles(j,i))
          end do
         end do
         Nrad = Nrad + Nregion(j)
        end do

        end if        
     
        close(13)

        end program XYZTS_generator  

!----------------------------------------------------------------------
!
!       Write out probe coordinates 
!
!----------------------------------------------------------------------
        subroutine write_probe(npr)

        use xyztsInp
        
        integer k, npr

        do k = 1, npr
           write(13,*) probe(k,1)-0.04, probe(k,2), probe(k,3)
           write(13,*) probe(k,1), probe(k,2), probe(k,3)
           write(13,*) probe(k,1)+0.02, probe(k,2), probe(k,3)
           write(13,*) probe(k,1)+0.04, probe(k,2), probe(k,3)
        enddo


        end subroutine 
