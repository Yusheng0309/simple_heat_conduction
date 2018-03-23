      program heatequ_2D

	use LinearAlgebra
	implicit none
c
c      data input
c      Lx =
c      Ly = 
c      dx = 
c      dy = 
c      Nx =
c      Ny =
c      total N
c      total t 
c      dt
c      output dt
c      T1~TN
c       111(i=0,j=10)  ---->     121(i=10,j=10)
c       -                        - 
c       -                        -
c       -                        -
c       1(i=0,j=0)     ---->     11(i=10,j=0)
c         
c        k = (j-1)*Nx + i
c        domain starting point (x0,y0)
c        x(k) = (i-1)*dx + x0
c        y(k) = (j-1)*dy + y0


c     dT/dt = D ( ddT/dx^2 + ddT/dy^2)
c     D = 
c     time derivative: forward finite difference
c     space derivative: Crank-Nicholson difference
c     a = D*dt/(dx^2)
c     b = D*dt/(dy^2)
c     (Nx-1)*(Ny-1) points have to be solved
 
c     boundary condition
c     T(y=Ly)=100
c     dT/dx(x=0,Lx , y=0)=0

c     initial condition
c     T(y=Ly,t)=100
c     other points = 0
       
	!matrix                     !N = (Nx-1)*(Ny-1)
	 integer:: N !inner points,total 81 points have to be solved
       integer:: p !total 121 grids including boundaries
	 real, allocatable :: A(:,:), B(:,:)      !p = (Nx+1)*(Ny+1) 
	 real, allocatable :: S(:), verb(:), c(:), z(:)
	 real, allocatable :: x(:), y(:), Tk(:)

      !numerical model parameters
       real:: Lx, Ly, Nx, Ny !Nx = Lx / dx = 10 (0~10, total 11)
	 real:: dx, dy, t, dt
	 real:: c_a, c_b
	 real:: D
	 integer:: outdt !output time
	 integer:: i, j, k, w, m !m=timestep
	 character supp*4,name*40
	 
      !domain length and width    
	 Lx=10.0
	 Ly=10.0
	!grid size and the node number
	 dx=0.25 !Lx=10
	 dy=0.25 !Ly=10
	 Nx=Lx/dx
	 Ny=Ly/dy
	 N = (Nx-1)*(Ny-1)
       p = (Nx+1)*(Ny+1) 
	 allocate( A(N,N), B(N,N) )     
	 allocate( S(N), verb(N), c(N), z(N) )
	 allocate( x(p), y(p), Tk(p) )

	!diffusion coefficient
	 D=1.0 ! t~sqrt(Ly/D)=10s !D = 0.01, 0.1 is OK

      !simulation time and output time
c	 t= 
	 dt=0.0625
c      m=t/dt         !timestep
       outdt=16       !...steps output once
	!heat equation discretization parameters
	 
	 c_a=D*dt/(dx**2)
	 c_b=D*dt/(dy**2)

      !initialization
	 do i=1,N
	   A(i,1:N)=0.0
	   B(i,1:N)=0.0
	   c(i)=0.0
       enddo


	!matrix A, coefficient matrix at next timestep 
	 do j = 1, Ny-1
	  do i = 1, Nx-1
	   k=(j-1)*(Nx-1)+i !(j-1)*(Nx-1)
	   A(k,k) = 2 + 2*c_a + 2*c_b
	   if( (i-1) > 0) then
	     A(k,k-1) = -c_a
	   else
	     A(k,k) = A(k,k) + (-c_a) !Neumann BC: left
	   endif
	   if( (i+1) < Nx) then ! <Nx
	     A(k,k+1) = -c_a
	   else
	     A(k,k) = A(k,k) + (-c_a) !Neumann BC: right
	   endif
	   if( (j-1) > 0) then
	     A(k,k-(Nx-1)) = -c_b
	   else
	     A(k,k) = A(k,k) + (-c_b) !Neumann BC: bottom
	   endif
	   if( (j+1) < Ny) then ! <Ny
	     A(k,k+(Nx-1)) = -c_b  !Dirichlet BC: top
	   endif 
        enddo
	 enddo  

	!matrix B, coefficient matrix at current timestep 
	 do j = 1, Ny-1
	  do i = 1, Nx-1
	   k=(j-1)*(Nx-1)+i !(j-1)*(Nx-1)
	   B(k,k) = 2 - 2*c_a - 2*c_b
	   if( (i-1) > 0) then
	     B(k,k-1) = c_a
	   else
	     B(k,k) = B(k,k) + (c_a) !Neumann BC: left
	   endif
	   if( (i+1) < Nx) then ! <Nx
	     B(k,k+1) = c_a
	   else
	     B(k,k) = B(k,k) + (c_a) !Neumann BC: right
	   endif
	   if( (j-1) > 0) then
	     B(k,k-(Nx-1)) = c_b
	   else
	     B(k,k) = B(k,k) + (c_b) !Neumann BC: bottom
	   endif
	   if( (j+1) < Ny) then ! <Ny
	     B(k,k+(Nx-1)) = c_b  !Dirichlet BC: top
	   endif 
        enddo
	 enddo  
           
      ! [A]S'=b=[B]S+c
	!vector c
	 do i = N-(Nx-1)+1, N !Dirichlet BC: top T=100 at y=Ly !N-(Nx-1)+1, N ,the last line
	  c(i) = 2*c_b*100
	 enddo

	! initial condition, S0
	 do i=1,N
	 S(i) = 0.0 !(T(i=1~(Nx-1) , j=1~(Ny-1))=0 , T(i,j=Ny)=100) 
	 enddo

      !coordinates of all grid points
	do j=0,Ny
	 do i=0,Nx
	  k=j*(Nx+1) + i +1 !No. 1~121 (1~(Nx+1)*(Ny+1))
	  x(k)= i*dx
	  y(k)= j*dy
	 enddo
	enddo
      
	
	!Solving loop-----------------------------------
       do m = 1,2000 ! t=10 dt

	  do i=1,N
	    verb(i)=0.0 
        enddo

	!b=[B]S+c
	 do i=1,N
	   do j=1,N
	     verb(i) =verb(i) + B(i,j)*S(j) 	     
	   enddo
	 enddo

	 do k=1,N
	  verb(k) = verb(k)+c(k)
	 enddo

	 call Gauss_Jordan(A, verb, S)

      !T:(Nx+1)*(Ny+1) points, S:(Nx-1)*(Ny-1) points
      !T(k) <- S(k)
	 do j=1,Ny-1
	  do i=1,Nx-1
	   k=j*(Nx+1) + i + 1 !No. 1~121 (1~(Nx+1)*(Ny+1)) 
	   w=(j-1)*(Nx-1) + i
	   Tk(k) = S(w)
	  enddo
	 enddo

	!T(k) for boundary points
	  do i=0,Nx
	  j=Ny !top boundary !No. 111~121
	  k=j*(Nx+1) + i + 1
        Tk(k)=100.
	  enddo
 
	  do j=1,Ny   !left up corner point is a repeat point
	  i = 0 !left boundary
	  k=j*(Nx+1) + i + 1
	  Tk(k)=Tk(k+1)
        enddo
	 
	  do j=1,Ny   !right up corner point is a repeat point
	  i = Nx !right boundary
	  k=j*(Nx+1) + i + 1
	  Tk(k)=Tk(k-1)
        enddo

        do i=0,Nx 
	  j=0 !bottom boundary !No. 1~11
	  k=j*(Nx+1) + i + 1
        Tk(k)=Tk(k+Nx+1)
	  enddo


      !output
	  if( mod(m,outdt) == 0. ) then !every outdt steps output once
        write(supp,'(i4.4)') m
	  name='PART_'//supp
	  write(*,*) name
        open(23,file=name,status='replace')
   
        call poute(22, m, x, y, Tk, p)

	  close(23)
	  endif

	 enddo !end solving loop----------------

	 stop
	end program

