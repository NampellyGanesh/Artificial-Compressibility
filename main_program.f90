program surface

use variables
use geometrical_properties 
use flux_calculation

implicit none
  !i face in vertical and j face in horizontal direction

 open(15,file="input.dat")

      read(15,*)
      read(15,*)
      read(15,*)file_name
      read(15,*)
      read(15,*)
      read(15,*)imax,jmax
      read(15,*)
      read(15,*)
      read(15,*)uinf,vinf,pinf,rho,c,lbc,rbc,tbc,bbc,istop,dt,restart
      read(15,*)
      read(15,*)
      read(15,*)CFL,irest,iprint,ipurge,beta

 close(15)
          
     
 allocate(darray(imax,jmax,2))
 allocate(i_norm(imax,2:jmax,2))
 allocate(j_norm(2:imax,jmax,2))
 allocate(cell_c(imax,jmax,2))

 allocate(iflux_sum(imax,2:jmax,3))
        
 allocate(jflux_sum(2:imax,jmax,3))
 allocate(sum_flux(2:imax,2:jmax,3))
 allocate(i_length(imax,2:jmax))
 allocate(j_length(2:imax,jmax))

 allocate(volume(2:imax,2:jmax))
 allocate(u(imax+1,jmax+1))
 allocate(v(imax+1,jmax+1))
 allocate(p(imax+1,jmax+1))
 allocate(delta_t(imax,jmax))
                                    
 print*,"allocated"

call i_normal()
call j_normal()
call volume_cell()
call cell_centre()          
call fluid_properties()



  deallocate(darray)
  deallocate(i_norm)
  deallocate(j_norm)
  deallocate(i_length)
  deallocate(j_length)
  deallocate(iflux_sum)
  deallocate(jflux_sum)
  deallocate(sum_flux)
!  deallocate(volume)
  deallocate(cell_c)
!  deallocate(delta_t)

 
print*,'deallocated'


end program surface
