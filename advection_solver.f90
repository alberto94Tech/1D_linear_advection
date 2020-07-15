program advection_solver
implicit none

integer :: i, N, u, numerical_scheme, test_case
real, dimension(:), allocatable :: x_coord, f0, f1, fa_0, fa_1, u_current_upwind, u_next_upwind, u_current_central, u_next_central, u_current_LW, u_next_LW, u_current_L, u_next_L, u_old, u_new, u_macc
real :: x_max, x_min, dx, t_current, t_final, dt, t, k, cfl

!Advection parameters and numerical option

open(1,file="input.par",status="old")
read(1,*) u !wave speed
read(1,*) x_min
read(1,*) x_max
read(1,*) N !number of grid points
read(1,*) cfl
read(1,*) numerical_scheme
read(1,*) test_case
 close(1)

!Array allocation

allocate (x_coord(N))
allocate (f0(N))
allocate (u_current_upwind(N))
allocate (u_next_upwind(N))
allocate (fa_0(N))
allocate (fa_1(N))
allocate (u_current_central(N))
allocate (u_next_central(N))
allocate (u_next_LW(N))
allocate (u_current_LW(N))
allocate (u_next_L(N))
allocate (u_current_L(N))
allocate (f1(N))
allocate (u_new(N))
allocate (u_old(N))
allocate (u_macc(N))

! grid generation
x_coord = 0.0
dx=(x_max-x_min)/(N-1) 
do i=1,N	
	if(i.eq.1)then
		x_coord(i)= x_min
	else
		x_coord(i) = x_min+dx*(i-1)	
	end if
end do

!time initialization

t_current = 0.0
t_final = 5.0

! time step
dt = cfl*dx/u

!selection of 
if(test_case.eq.0)then

	!initial condition for f0
	f0=0.0
	do i=1, N
		f0(i)=0.5*(sign(1.0, x_coord(i))+1)
	end do

	! analytical solution fa_0
	do i=1, N

		fa_0(i)=0.5*(sign(1.0, x_coord(i)-1.5*t)+1)
	
	end do

	! boundary conditions for f0

	u_next_upwind(1) = 0.0
	u_next_upwind(N) = 1.0

	u_next_central(1) = 0.0
	u_next_central(N) = 1.0

	u_next_L(1) = 0.0
	u_next_L(N) = 1.0

	u_next_LW(1) = 0.0
	u_next_LW(N) = 1.0	

	u_new(1) = 0.0
	u_new(N) = 1.0

	u_macc(1) = 0.0
	u_macc(N) = 1.0	

	! initialization of the numerical scheme 
	u_current_upwind = f0
	u_current_central = f0
	u_current_LW = f0
	u_current_L = f0
	u_old = f0

else if(test_case.eq.1)then

	! Initial condition f1

	do i=1, N
		f1(i)=0.5*(exp(-x_coord(i)*x_coord(i)))
	end do

	! analytical solution fa_1

	do i=1, N

		fa_1(i)=0.5*(exp(-(x_coord(i)-1.5*t)**2))
	end do

	! boundary conditions for f1

	u_next_upwind(1) = 0.0
	u_next_upwind(N) = 0.0

	u_next_central(1) = 0.0
	u_next_central(N) = 0.0

	u_next_L(1) = 0.0
	u_next_L(N) = 0.0

	u_next_LW(1) = 0.0
	u_next_LW(N) = 0.0
	
	u_new(1) = 0.0
	u_new(N) = 0.0

	u_macc(1) = 0.0
	u_macc(N) = 0.0	

	! initialization of the numerical scheme 

	u_current_upwind = f1
	u_current_central = f1
	u_current_LW = f1
	u_current_L = f1
	u_old = f1


else
	print*,'Error'
end if

! numerical solutions
do while(t_current.lt.t_final)
	
	do i=1, N-1
		u_next_upwind(i) = u_current_upwind(i) - (u*dt/dx)*(u_current_upwind(i) - u_current_upwind(i-1))
		u_next_central(i) = u_current_central(i)-(u*dt/dx)*(u_current_central(i+1)- u_current_central(i-1))/2.0
		u_next_L(i) = 0.5*(u_current_L(i+1)+u_current_L(i-1))-(u*dt/2*dx)*(u_current_L(i+1)-u_current_L(i-1))
		u_next_LW(i) = u_current_LW(i)-((u*dt)/(2.0*dx))*(u_current_LW(i+1)-u_current_LW(i-1))+0.5*u**2*dt**2*(1/(dx**2))*(u_current_LW(i+1)-2*u_current_LW(i)+u_current_LW(i-1))
		u_macc(i)=u_old(i)-(u*dt/dx)*(u_old(i+1)-u_old(i))
		u_new(i)=0.5*(u_old(i)+u_macc(i))-(u*dt/dx)*(u_macc(i)-u_macc(i-1))
		
	end do

	u_current_upwind = u_next_upwind
	u_current_central = u_next_central
	u_current_LW = u_next_LW
	u_current_L = u_next_L
	u_old=u_new

	t_current=t_current+dt

end do

!write results on a file

if(test_case.eq.0)then

	open (unit=2,file="f0_1.dat",status="replace",action="write")
 		write(2,*) 'TITLE="Numerical and Analytical Solution for f_0"'
  		write(2,*) 'VARIABLES= "X" "ANALYTICAL" "u_UPWIND" "u_CENTRAL"'
 		write(2,*) "ZONE I=",N+1," F=POINT"

		do i=1,N
			write(2,*) x_coord(i),fa_0(i), u_next_upwind(i), u_next_central(i) 
		end do
	close (2)


	open(unit=3,file="f0_2.dat",status="replace",action="write")
		write(3,*) 'TITLE="Numerical and Analytical Solution for f_0"'
  		write(3,*) 'VARIABLES= "X" "u_LAX-FRIEDRICHS" "u_LAX-WENDROFF" "u_MACCORMACK" '
 		write(3,*) "ZONE I=",N+1," F=POINT"

		do i=1, N
			write(3,*)  x_coord(i), u_next_L(i), u_next_LW(i),u_new(i) 			
		end do

	close(3)

else if(test_case.eq.1)then

open (unit=4,file="f1_1.dat",status="replace",action="write")
 		write(4,*) 'TITLE="Numerical and Analytical Solution for f_1"'
  		write(4,*) 'VARIABLES= "X" "ANALYTICAL" "u_UPWIND" "u_CENTRAL"'
 		write(4,*) "ZONE I=",N+1," F=POINT"

		do i=1,N

			write(4,*) x_coord(i),fa_1(i), u_next_upwind(i), u_next_central(i) 
		end do
	close (4)

	open(unit=5,file="f1_2.dat",status="replace",action="write")
		write(5,*) 'TITLE="Numerical and Analytical Solution for f_1"'
  		write(5,*) 'VARIABLES= "X" "u_LAX-FRIEDRICHS" "u_LAX-WENDROFF" "u_MACCORMACK" '
 		write(5,*) "ZONE I=",N+1," F=POINT"

		do i=1, N
			write(5,*) x_coord(i), u_next_L(i), u_next_LW(i), u_new(i) 			
		end do
	close(5)

end if
end program advection_solver
