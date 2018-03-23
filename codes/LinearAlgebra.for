      module LInearAlgebra
	 implicit none
	contains
	
	!Gauss_Jordan elimination 
	subroutine Gauss_Jordan(A, S, ANS)
	 
	 implicit none
	 real :: A(:,:)
	 real :: S(:)
	 real :: ANS(:)
	 real, allocatable :: B(:,:)
	 integer :: i, N
	 N = size(A,1)
	 allocate( B(N,N) )
	 !retaining original matrix A and array S
	 B = A
	 ANS = S
	 !Diagonalizing B to a diagonal matrix
	 call Upper(B, ANS, N) ! simplifying B to a upper triangular matrix
	 call Lower(B, ANS, N) ! Then simplifying B to a lower triangular matrix

	 !solving
	 forall(i=1:N)
	  ANS(i) = ANS(i)/B(i,i)
	 endforall
	 return
	end subroutine Gauss_Jordan
	 
	!output linear algebric equations 
	subroutine output(M,S)
	 implicit none
	 real :: M(:,:), S(:)
	 integer :: N, i, j
	 N = size(M,1)
	 !advance="no" is to prevent from line breaks, and continue the line kept in the same raw
	 do i = 1, N
	   write(*,"(1x,f5.2,a1)", advance="NO") M(i,1), 'A' !element(1,1) has to be positive
	   do j = 2, N
	      if( M(i,j) < 0 ) then
	         write(*,"('-',f5.2,a1)", advance="NO") -M(i,j),char(64+j)
            else 
	         write(*,"('+',f5.2,a1)", advance="NO") M(i,j),char(64+j)
	      endif
         enddo
	   write(*,"('=',f8.4)") S(i)
	  enddo
	  return
	end subroutine output

      !output matrix mn
	subroutine output_matrix(matrix,b)
	  implicit none
	  integer :: m, n
	  real :: matrix(:,:), b(:)
	  integer :: i
	  character(len=20) :: for='(??(1x,f8.4))'
	  m = size(matrix, 1)
	  n = size(matrix, 2)
        
	  open(10, file="matrix.txt")
        write( FOR(2:3), '(I2)' ) N
	  do i = 1, N
	    write(10,FMT=FOR ) matrix(i,:)
	  enddo
	  close(10)

c	  write(*,*) 'vector b:'
c	 do i = 1, N
c	  write(*,"(1x,f10.4)") b(i)
c	 enddo
	  return
	end subroutine output_matrix

	     
      ! Upper
	subroutine Upper(M, S, N)
	 implicit none
	 integer :: N
	 real :: M(N,N)
	 real :: S(N)
	 integer :: I, J
	 real :: EU
	 do I = 1, N-1
	   do J = I+1, N
	   EU = M(J,I)/M(I,I)
	   M(J,I:N) = M(J,I:N)-M(I,I:N)*EU
	   S(J) = S(J)-S(I)*EU
	   enddo
	 enddo
	 return
	end subroutine Upper
	
	! Lower
	subroutine Lower(M, S, N)
	 implicit none
	 integer :: N
	 real :: M(N,N)
	 real :: S(N)
	 integer :: I, J
	 real :: EL
	 do I = N, 2, -1
	  do J = I-1, 1, -1
	   EL = M(J,I)/M(I,I)
	   M(J,1:N) = M(J,1:N)-M(I,1:N)*EL   
	   S(J) = S(J)-S(I)*EL
        enddo
	 enddo
	 return 
	end subroutine
	end module