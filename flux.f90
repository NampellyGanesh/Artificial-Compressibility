module flux_calculation
   use variables
   use geometrical_properties
contains
        subroutine fluid_properties
           implicit none   
!           beta=scale_fac*uinf
         !intialize the flow field
         if(restart .eq. 0)then
              do j=2,jmax 
                 do i=2,imax
                    u(i,j)=uinf
                    v(i,j)=vinf
                    p(i,j)=pinf
                 end do
              end do
              residue=(/0.0,0.0,0.0/) 
          !restart the file
          else
               write(filename_results,120)ipurge
               open(40,file=filename_results)
               print*,'restart' 
               read(40,*)
               read(40,*)
               do j=2,jmax
                   do i=2,imax
                      read(40,*)dummy1,dummy2,u(i,j),v(i,j),p(i,j)
                    end do
               end do
               close(40)
           end if
           no_of_cells=(imax-1)*(jmax-1)  
 
           do iteration=1,istop    

          ! Enforce Boundary conditions on ghost cells
          !    1      =     inlet
          !    2      =     outlet
          !    3      =     slip wall
          !    4      =     No slip wall
   
          if(lbc.eq.1) then
                      u(1,2:jmax)=uinf
                      v(1,2:jmax)=vinf
                      p(1,2:jmax)=p(2,2:jmax)

           else if (lbc.eq.2) then
                      u(1,2:jmax)=u(2,2:jmax)
                      v(1,2:jmax)=v(2,2:jmax)
                      p(1,2:jmax)=pinf

            else if (lbc.eq.3) then
                      u(1,2:jmax)=u(2,2:jmax)-2*(u(2,2:jmax)*i_norm(1,2:jmax,1)&
                                     +v(2,2:jmax)*i_norm(1,2:jmax,2))*i_norm(1,2:jmax,1)
                      v(1,2:jmax)=v(2,2:jmax)-2*(u(2,2:jmax)*i_norm(1,2:jmax,1)&
                                     +v(2,2:jmax)*i_norm(1,2:jmax,2))*i_norm(1,2:jmax,2)
                      p(1,2:jmax)=p(2,2:jmax)
            
            else if (lbc.eq.4) then
                      u(1,2:jmax)=-u(2,2:jmax)  
                      v(1,2:jmax)=-v(2,2:jmax)
                      p(1,2:jmax)=p(2,2:jmax)
           
            end if
        !Right boundary
       
            if(rbc.eq.1) then
                      u(imax+1,2:jmax)=uinf
                      v(imax+1,2:jmax)=vinf
                      p(imax+1,2:jmax)=p(imax,2:jmax)
                                                                                              
           else if (rbc.eq.2) then
                      u(imax+1,2:jmax)=u(imax,2:jmax)
                      v(imax+1,2:jmax)=v(imax,2:jmax)
                      p(imax+1,2:jmax)=pinf
                                                                                              
           else if (rbc.eq.3) then
                      u(imax+1,2:jmax)=u(imax,2:jmax)-2*(u(imax,2:jmax)*i_norm(imax,2:jmax,1)&
                                       +v(imax,2:jmax)*i_norm(imax,2:jmax,2))*i_norm(imax,2:jmax,1)
                      v(imax+1,2:jmax)=v(imax,2:jmax)-2*(u(imax,2:jmax)*i_norm(imax,2:jmax,1)&
                                       +v(imax,2:jmax)*i_norm(imax,2:jmax,2))*i_norm(imax,2:jmax,2)
                      p(imax+1,2:jmax)=p(imax,2:jmax)
                                                                                              
           else if (rbc.eq.4) then
                     u(imax+1,2:jmax)=-u(imax,2:jmax)  
                     v(imax+1,2:jmax)=-v(imax,2:jmax)
                     p(imax+1,2:jmax)=p(imax,2:jmax)
           end if 

             !top boundary
           if(tbc.eq.1) then                                                                                                                               
                     u(2:imax,jmax+1)=uinf
                     v(2:imax,jmax+1)=vinf
                     p(2:imax,jmax+1)=p(2:imax,jmax)
                                                                                               
            else if (tbc.eq.2) then
                     u(2:imax,jmax+1)=u(2:imax,jmax)
                     v(2:imax,jmax+1)=v(2:imax,jmax)
                     p(2:imax,jmax+1)=pinf
            else if (tbc.eq.3) then
                                                                                               
                     u(2:imax,jmax+1)=u(2:imax,jmax)-2*(u(2:imax,jmax)*j_norm(2:imax,jmax,1)&
                                      +v(2:imax,jmax)*j_norm(2:imax,jmax,2))*j_norm(2:imax,jmax,1)
                     v(2:imax,jmax+1)=v(2:imax,jmax)-2*(u(2:imax,jmax)*j_norm(2:imax,jmax,1)&
                                      +v(2:imax,jmax)*j_norm(2:imax,jmax,2))*j_norm(2:imax,jmax,2)
                     p(2:imax,jmax+1)=p(2:imax,jmax)                                                                                                          
            else if (tbc.eq.4) then
                                                                                               
                     u(2:imax,jmax+1)=-u(2:imax,jmax)  
                     v(2:imax,jmax+1)=-v(2:imax,jmax)
                     p(2:imax,jmax+1)=p(2:imax,jmax)
            end if 
            
           !Bottom Boundary

           if(bbc.eq.1) then                                                                                                                               
                     u(2:imax,1)=uinf
                     v(2:imax,1)=vinf
                     p(2:imax,1)=p(2:imax,2)
                                                                                                    
           else if (bbc.eq.2) then
                     u(2:imax,1)=u(2:imax,2)
                     v(2:imax,1)=v(2:imax,2)
                     p(2:imax,1)=pinf
                                                                                                    
           else if (bbc.eq.3) then
                     u(2:imax,1)=u(2:imax,2)-2*(u(2:imax,2)*j_norm(2:imax,1,1)&
                                   +v(2:imax,2)*j_norm(2:imax,1,2))*j_norm(2:imax,1,1)
                     v(2:imax,1)=v(2:imax,2)-2*(u(2:imax,2)*j_norm(2:imax,1,1)&
                                  +v(2:imax,2)*j_norm(2:imax,1,2))*j_norm(2:imax,1,2)
                     p(2:imax,1)=p(2:imax,2)                                                                                                          
           else if (bbc.eq.4) then
                     u(2:imax,1)=-u(2:imax,2)  
                     v(2:imax,1)=-v(2:imax,2)
                     p(2:imax,1)=p(2:imax,2)
           end if           
                     
          !calculate flux
          !I direction flux in horizontal

          do j=2,jmax
             do i=1,imax   
                     l_velocity(:)=(/1.0,u(i,j),v(i,j)/)
                     r_velocity(:)=(/1.0,u(i+1,j),v(i+1,j)/)
                     inormal1(:)=(/0.0,i_norm(i,j,1),i_norm(i,j,2)/)
                     pressure_gradient=p(i,j)-p(i+1,j)

                     x=0.5*((u(i,j)+u(i+1,j))*i_norm(i,j,1)+(v(i,j)+v(i+1,j))*i_norm(i,j,2))
                    lambda=(abs(x)+sqrt(x**2+4*beta**2))*0.5
                    
                    fluxl=max(0.0,(u(i,j)*i_norm(i,j,1)+v(i,j)*i_norm(i,j,2)))
                    fluxr=min(0.0,(u(i+1,j)*i_norm(i,j,1)+v(i+1,j)*i_norm(i,j,2)))
                    ifluxl(:) = rho*(fluxl+((c*pressure_gradient)/(2*rho*lambda)))*l_velocity(:)
                    ifluxr(:) = rho*(fluxr+((c*pressure_gradient)/(2*rho*lambda)))*r_velocity(:)

            !calculating pressure flux
                  
                    pressure_interface=0.5*(p(i,j)+p(i+1,j))
                    pressure_flux(:)=pressure_interface*inormal1(:)
                    iflux_sum(i,j,:)=(ifluxl(:)+pressure_flux(:)+ifluxr(:))*i_length(i,j)
            end do

         end do

                 !J faces FLUX IN VERTICAL 

         do j=1,jmax
            do i=2,imax 

                    r_velocity(:)=(/1.0,u(i,j+1),v(i,j+1)/)
                    l_velocity(:)=(/1.0,u(i,j),v(i,j)/) 
                    jnormal1(:)=(/0.0,j_norm(i,j,1),j_norm(i,j,2)/)
                    pressure_gradient=(p(i,j)-p(i,j+1))

                    x=0.5*((u(i,j)+u(i,j+1))*j_norm(i,j,1)+(v(i,j)+v(i,j+1))*j_norm(i,j,2))
                    lambda=(abs(x)+sqrt(x**2+4*beta**2))*0.5
                    fluxb=max(0.0,(u(i,j)*j_norm(i,j,1)+v(i,j)*j_norm(i,j,2)))
                    fluxt=min(0.0,(u(i,j+1)*j_norm(i,j,1)+v(i,j+1)*j_norm(i,j,2)))

                    jfluxb(:)=rho*(fluxb+((c*pressure_gradient)/(2*lambda*rho)))*l_velocity(:)                      
                    jfluxt(:)=rho*(fluxt+((c*pressure_gradient)/(2*lambda*rho)))*r_velocity(:)
                
                    !pressure flux
                   
                    pressure_interface=0.5*(p(i,j+1)+p(i,j))
                    pressure_flux(:)=pressure_interface*jnormal1(:)
  !                  if(j.eq.1 .or. j.eq.jmax) then
  !                     jfluxb(:) = 0.0 
  !                     jfluxt(:) = 0.0
  !                  end if 

                    jflux_sum(i,j,:)=(jfluxb(:)+pressure_flux(:)+jfluxt(:))*j_length(i,j)
           end do 

        end do

        !local time stepping
        do j=2,jmax
           do i=2,imax
                    !finding averages
                    u_int1=(u(i,j)+u(i+1,j))*0.5
                    u_int2=(u(i,j)+u(i-1,j))*0.5
                    u_int3=(u(i,j)+u(i,j+1))*0.5
                    u_int4=(u(i,j)+u(i,j-1))*0.5

                    v_int1=(v(i,j)+v(i+1,j))*0.5    
                    v_int2=(v(i,j)+v(i-1,j))*0.5
                    v_int3=(v(i,j)+v(i,j+1))*0.5
                    v_int4=(v(i,j)+v(i,j-1))*0.5

                    u1=u_int1*i_norm(i,j,1)+v_int1*i_norm(i,j,2)
                    u2=u_int2*i_norm(i-1,j,1)+v_int2*i_norm(i-1,j,2)
                    u3=u_int3*j_norm(i,j,1)+v_int3*j_norm(i,j,2)
                    u4=u_int4*j_norm(i,j-1,1)+v_int4*j_norm(i,j-1,2)
           
                   !right face time step
                 
                    lamda1=0.5*(abs(u1)+sqrt(u1**2+4*beta**2))
                             
                   !left face time step
           
                    lamda2=0.5*(abs(u2)+sqrt(u2**2+4*beta**2))
             
                   !top face time step
           
                    lamda3=0.5*(abs(u3)+sqrt(u3**2+4*beta**2)) 
           
                    !bottom face time step
              
                    lamda4=0.5*(abs(u4)+sqrt(u4**2+4*beta**2)) 
            
                    delta_t(i,j)=(0.5*(i_length(i,j)*lamda1+i_length(i-1,j)*lamda2+&
                                 j_length(i,j)*lamda3+j_length(i,j-1)*lamda4))/CFL
           end do
        end do              

        !calculate sum of fluxes at cell centre
        do j=2,jmax
           do i=2,imax  
                    sum_flux(i,j,:)=iflux_sum(i,j,:)-iflux_sum(i-1,j,:)+jflux_sum(i,j,:)-jflux_sum(i,j-1,:)
                    if(isnan(sum_flux(i,j,1)))then
                      print*,'NAN is encountered'
                      stop
                    end if
                    Q1(:)=(/p(i,j),u(i,j),v(i,j)/)
                    Q1(:) = Q1(:)-(/beta**2*sum_flux(i,j,1),&
                                       -(u(i,j)/rho)*sum_flux(i,j,1)+sum_flux(i,j,2)/rho ,&
                                        -(v(i,j)/rho)*sum_flux(i,j,1)+(sum_flux(i,j,3)/rho)/)/delta_t(i,j)
                    p(i,j)=Q1(1)
                    u(i,j)=Q1(2)
                    v(i,j)=Q1(3)
                           
                    residue(:)=residue(:)+abs(sum_flux(i,j,:))
            end do
        end do

        residue(:) = residue(:)/no_of_cells
        open(2,file='results/residual.dat')
        print*,residue(:),iteration,maxval(delta_t(:,:))
        write(2,100)iteration,residue
        100 format(1i10,3e30.15)
        120 format('results/ycuttcell',i1,'.dat')                                           
        if(mod(iteration,irest) .eq. 0.0)then
              write(filename_results,120)h
              open(16,file=filename_results)
              write(16,*)'Variables=   x, y, u , v, p,  rho'
              write(16,*)'zone i=', imax-1,'    j='   ,  jmax-1
              do j=2,jmax                                                                 
                                                                      
                 do i=2,imax
                                                                           
                    write(16,*)cell_c(i,j,1:2),u(i,j),v(i,j),p(i,j)
                                                                           
                  end do
               end do
               h=h+1
               print*,'continuity,        x-momentum      y-momentum      iteration'
               close(16)
               if(h .gt. iprint)then
                  h=1
               end if
               else 
                  continue
         end if 
         if (maxval(residue).lt.1.0e-8)then
               print*, "solution is converged"
               goto 110
         else 
               continue                                        
         end if
       end do
       close(2)
       110  write(filename_results,120)h
       open(16,file=filename_results)
       write(16,*)'Variables=   x, y, u , v, p'
       write(16,*)'zone i=', imax-1,'    j='   ,  jmax-1

       do j=2,jmax
          do i=2,imax
             write(16,*)cell_c(i,j,1:2),u(i,j),v(i,j),p(i,j)
          end do
       end do
       close(16)

end subroutine fluid_properties

end module flux_calculation




