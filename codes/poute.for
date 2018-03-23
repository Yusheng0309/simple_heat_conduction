      subroutine poute(nf, m, x, y, Tk, p)
      use LinearAlgebra
      
      integer nf, i, p, m
	real :: x(p), y(p), Tk(p)
      !m: timesteps, p: total grid points
      do i=1,p
         write(nf+1,100) x(i), y(i), Tk(i), i
      enddo

100   format(3e16.8,i8)

c     -------------------------------------

       write(*,*) 't=',m,' steps','  ok'


      return 
      end

