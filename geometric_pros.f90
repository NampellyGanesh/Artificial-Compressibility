module geometrical_properties
use variables

   implicit none
   ! i face in vertical and j  face in horizontal direction
   !length of i faces

contains 
     subroutine i_normal()

   ! print*,"Enter the file name"
   !  read*,mesh_file

   open(10,file=file_name)
   n=imax*jmax
   do j=1,jmax
      do i=1,imax
         read(10,*)darray(i,j,1:2)
      end do
   end do
   close(10)

   print*,'geometric_variables'     
   do j=2,jmax
      do i=1,imax
         i_length(i,j)=sqrt((darray(i,j,1)-darray(i,j-1,1))**2+(darray(i,j,2)-darray(i,j-1,2))**2)  !length of I face

         !normals of i face
         i_res(1:2)=darray(i,j,1:2)-darray(i,j-1,1:2)
         i_norm(i,j,1)=k(2)*0.0+k(3)*i_res(2)
         i_norm(i,j,2)=-(k(1)*0.0+k(3)*i_res(1))
         magnitude=sqrt(i_norm(i,j,1)**2+i_norm(i,j,2)**2)
         i_norm(i,j,:)=i_norm(i,j,:)/magnitude
      end do 
   end do
   !writing into file 
   open(11,file='results/i_length.dat')
   write(11,*)'length, xnormal, ynormal, x,y'

   do j=2,jmax
      do i=1,imax
         write(11,'(7e30.15)')i_length(i,j),i_norm(i,j,1)/50.0,i_norm(i,j,2)/50.0,darray(i,j,1),darray(i,j,2)
      end do
   end do
          
   close(11)
      
   end subroutine i_normal

   subroutine j_normal()
 
     do j=1,jmax
        do i=2,imax
           j_length(i,j)=sqrt((darray(i,j,1)-darray(i-1,j,1))**2+(darray(i,j,2)-darray(i-1,j,2))**2) !J Face length
 
           !normals of j face
           j_res(1:2)=darray(i,j,1:2)-darray(i-1,j,1:2)
           j_norm(i,j,1)=k(2)*0.0-k(3)*j_res(2)
           j_norm(i,j,2)=-(k(1)*0.0-k(3)*j_res(1))
           magnitude=sqrt(j_norm(i,j,1)**2+j_norm(i,j,2)**2)
           j_norm(i,j,:)=j_norm(i,j,:)/magnitude
        end do
     end do
     !writing into file
     open(12,file='results/j_length.dat')
     write(12,*)'length,xnormal, ynormal,x,y'
 
     do j=1,jmax
        do i=2,imax
            write(12,'(7e30.15)')j_length(i,j),j_norm(i,j,1)/50.0,j_norm(i,j,2)/50.0,darray(i,j,1),darray(i,j,2) 
        end do
     end do
     close(12)

    end subroutine j_normal

    !!calculate volume of each cell

    subroutine volume_cell()

     do j=2,jmax
        do i=2,imax
           diag1(1:2)=darray(i,j,1:2)-darray(i-1,j-1,1:2)
           diag2(1:2)=darray(i-1,j,1:2)-darray(i,j-1,1:2)
           cross_diag=diag1(1)*diag2(2)-diag1(2)*diag2(1)
           volume(i,j)=abs(cross_diag)/2.0
        end do
     end do
        
     open(13,file='results/Volume_cell.dat')
     write(13,*)'i,j,volume'
     do j=2,jmax
        do i=2,imax
           write(13,*)i,j,volume(i,j)
        end do
     end do

     close(13)

    end subroutine volume_cell
     
!Calculate cell_Centres

    subroutine cell_centre()
   
    do j=2,jmax
       do i=2,imax
          cell_c(i,j,:)=(darray(i,j,:)+darray(i,j-1,:)+darray(i-1,j-1,:)+darray(i-1,j,:))/4.0
       end do
    end do

    open(14,file='results/cell_centre.dat')
    write(14,*)'i,j,cell_centre'
    do j=2,jmax
       do i=2,imax
           write(14,*),i,j,cell_c(i,j,1),cell_c(i,j,2)
       end do
    end do

   close(14)

   end subroutine cell_centre

end module geometrical_properties





