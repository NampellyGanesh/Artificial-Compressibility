module variables

       implicit none
       
       real :: cross_diag,magnitude,uinf,vinf,pinf,rho,scale_fac,beta,fluxl,fluxr,x,lambda,pressure_interface,dt,fluxb,fluxt
       real    ::     dummy1,dummy2,pressure_gradient,u1,u2,u3,u4,lamda1,lamda2,lamda3,lamda4,CFL
       real     ::    u_int1,u_int2,u_int3,u_int4,v_int1,v_int2,v_int3,v_int4,c
       integer :: i,j,imax,jmax,n,lbc,rbc,tbc,bbc,no_of_Cells,istop,iteration,restart,irest,h=1,iprint,ipurge
       character(len=40) :: file_name,filename_results
       real, dimension(:,:,:), allocatable    ::  darray,i_norm,j_norm,cell_c,iflux_sum,jflux_sum,sum_flux
       real, dimension(:,:),   allocatable    ::  i_length,j_length,volume,u,v,p,delta_t
       real, dimension(2)                     ::  diag1,diag2,j_res,i_res
       real, dimension(3)                     ::  pressure_flux,r_velocity,l_velocity,b_velocity,inormal1,jnormal1,&
                                                  ifluxl,ifluxr,jfluxb,jfluxt, k= [0,0,1] 
       real, dimension(3)                     ::  Q1,residue
       
contains
        subroutine allocation() 
          allocate(darray(imax,jmax,2))
          allocate(i_norm(imax,2:jmax,2))
          allocate(j_norm(2:imax,jmax,2))
          allocate(cell_c(imax,jmax,2))
          
          allocate(iflux_sum(imax,2:jmax,3))
                  
          allocate(jflux_sum(2:imax,jmax,3))
          allocate(sum_flux(2:imax,2:jmax,3))
          allocate(i_length(imax,2:jmax))
          allocate(j_length(2:imax,jmax))
          
          allocate(volume(imax,jmax))
          allocate(u(imax,jmax))
          allocate(v(imax,jmax))
          allocate(p(imax,jmax))
          allocate(delta_t(imax,jmax))

          print*,"allocated"
        end subroutine allocation


end module variables
