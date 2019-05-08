module task
use :: mpi
implicit none
contains

subroutine GetMaxCoordinates(A,x1,y1,x2,y2)
   real(8),dimension(:,:),intent(in) :: A
   integer(4),intent(out):: x1,x2,y1,y2
   integer(4) :: n, l, r, up, down, m, emp
   real(8), allocatable :: current_column(:), B(:,:)
   real(8) :: current_sum, max_sum
   logical :: transpos
   integer(4) :: mpiErr,mpiSize,mpiRank
   integer(4),dimension(MPI_STATUS_SIZE) :: status
   real(8),allocatable,dimension(:) :: Max_value_submatrixA
   real(8) numRank_with_maxsubA

   call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
   call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
   
   m = size(A,1) 
   n = size(A,2) 
   transpos=.false.
   
   if(m<n) then
     transpos=.true.
     B=transpose(A)
     m=size(B,1)
     n=size(B,2)
   else
     B=A
   endif

allocate(current_column(m))

max_sum=B(1,1)
x1=1
x2=1
y1=1
y2=2

   do l=mpiRank+1,n,mpiSize
     current_column=B(:,l)
     do r=l,n
         if(r>l) then
           current_column=current_column+B(:,r)  
         endif
      
     call GetMaxInArray(current_column,current_sum,up,down)
          
          if(current_sum>max_sum) then
            max_sum=current_sum
            x1=up
            x2=down
            y1=l
            y2=r
          endif
     enddo
   enddo

allocate(Max_value_submatrixA(0:mpiSize-1))

call mpi_gather(max_sum, 1, MPI_REAL8, Max_value_submatrixA(0:mpiSize-1),mpiSize , MPI_REAL8, 0,MPI_COMM_WORLD,mpiErr)

 if(mpiRank==0) then
   numRank_with_maxsubA=maxloc(Max_value_submatrixA(0:mpiSize-1),1)
   endif



deallocate(current_column)
 
  if(transpos) then
     emp=x1
     x1=y1
     y1=emp
     emp=y2
     y2=x2
     x2=emp
  endif
end subroutine


subroutine GetMaxInArray(c,sum,up,down)
   real(8), intent(in), dimension(:) :: c 
   integer(4), intent(out) :: up, down
   real(8), intent(out) :: sum
   real(8) :: cur_sum
   integer(4) :: minus_pos, i

   sum=c(1)
   up=1
   down=1
   cur_sum=0
   minus_pos=0

   do i=1,size(c)
      cur_sum=cur_sum+c(i)
       if(cur_sum>sum) then
          sum=cur_sum
          up=minus_pos+1
          down=i
       endif
       
       if(cur_sum<0) then
          cur_sum=0
          minus_pos=i
       endif
    enddo
end subroutine GetMaxInArray
end module task
         





















