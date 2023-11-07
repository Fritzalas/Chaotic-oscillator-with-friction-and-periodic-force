program rk_solve
  use constant_module                                 !In this module we have the values of PI and 2*PI and also for omega_0,omega,gamma,a_0,omega_02,omega2,Energy
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8),allocatable       :: T(:),X1(:),X2(:)       !We have the arrays of our problem. T is the array of time,X1 is the array for space and X2 is the array for velocity
  real(8)                   :: Ti,Tf                  !Initial time and final time
  real(8)                   :: X10,X20                !Initial values for space and velocity
  integer                   :: Nt                     !The total steps from ti to tf
  integer                   :: i,j                    !Useful for do loops
  real(8)                   :: t_trans                !The transient time that we are going to get and see our data after that
  integer                   :: index                  !A useful integer that shows that the transient behaviour has finished
  integer                   :: i_trans                !An integer that multiplies pi to give the t_trans
  integer                   :: i_tf                   !An integer that multiplies pi to give the tf
  integer                   :: ktrans                 !An integer that we are going to use for not transient times
  integer                   :: kmax                   !The maximum k that we are going to reach at our final time
  real(8)                   :: dummy_x1,dummy_x2      !The dummy indexes are going to hold the starting values of our problem
  real(8)                   :: dummy_t                !The dummy index that is going to hold the starting time every time for our problem
  real(8),dimension(50)     :: Q,Q1,Q2                !Dummy vectors for position,angular position and angular velocity
  real(8)                   :: help
  !We don't need them to declare because of the constant_module
  !real(8)                   :: omega_0,omega          !Angular velocity for the mass and also for the extra periodic force
  !real(8)                   :: omega_02,omega2        !The square of the above angular velocities
  !real(8)                   :: gamma,a_0              !The reduction variable and also the length of the periodic force
  !real(8)                   :: Energy                 !The total energy of our system
  !---------------------------------------------------------------------------------
  !User Interface:
  print *,'# Runge-Kutta Method for FDP Integration'
  print *,'# Enter omega_0, omega, gamma, A:'
  read  *,         omega_0, omega, gamma, a_0
  print *,'# omega_0= ',omega_0,' omega= ', omega
  print *,'# gamma=   ',gamma,  ' A=     ',a_0
  print *, "#Enter the number of spaces in [ti,tf]: "
  read  *, Nt
  print *, "#Enter the initial time ti and the index of the final time tf: "
  read  *, Ti,i_tf
  print *, "#Enter the initial values for space and velocity: "
  read  *, X10,X20
  print *, "#We have Nt,Ti,i_Tf,X10,X20 equal to: ",Nt,Ti,i_tf,X10,X20
  print *, "#Give to the program the index for the transient time t_trans:"
  read  *, i_trans
  print *, "#i_trans=",i_trans
  !----------------------------------------------------------------------------------
  !Array Allocation:
  ALLOCATE(X1(Nt))
  ALLOCATE(X2(Nt))
  ALLOCATE(T(Nt))
  !----------------------------------------------------------------------------------
  !Fix the times tf and t_trans with pi:
  Tf=i_tf*PI
  t_trans=i_trans*PI
  !----------------------------------------------------------------------------------
  !Mistake in given data:
  if (Tf<Ti) stop 'The final time is less than the initial time! Try Again!!!'
  !----------------------------------------------------------------------------------
  !Calculations:
  omega_02 = omega_0*omega_0
  omega2   = omega  *omega            !Calculate omega square for next calculations
  ktrans=i_trans
  kmax=i_tf                           !We give the right values to our problem
  call RK(T,X1,X2,Ti,Tf,X10,X20,Nt)   !Useful subroutine to make the calculations
  !----------------------------------------------------------------------------------
  !Transient behaviour:
  if (t_trans>tf) stop 'The transient time is bigger than the final time! Try Again!!!'
  do i=1,Nt
     if (T(i)>=t_trans) then
        index=i
        dummy_x1=X1(index)
        dummy_x2=X2(index)
        exit   !We found out our position of the after transient behaviour
     end if
     !-------------------------------------------------------------------------------
  end do
  dummy_t=i_trans*PI
  !----------------------------------------------------------------------------------
  !Output file:
  !Calculation of time,angular position,angular velocity and energy:
  open(unit=11,file='logistic_map.dat')
  do i=ktrans,kmax-1
     dummy_t=i*PI
     do j=index,Nt
        if (T(j)>=dummy_t) then
           write(11,*) a_0,T(j),X1(j),X2(j)
           exit
        end if
     end do
  end do
  close(11)
  !----------------------------------------------------------------------------------
end program rk_solve
!====================================================================================
real(8) function f1(t,x1,x2)
  use constant_module
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8)     :: t,x1,x2 !Given time,space and velocity
  !----------------------------------------------------------------------------------
  !Calculations:
  f1=x2     !Because dx/dt=V=x2 from our problem
end function f1
!====================================================================================
real(8) function f2(t,x1,x2)
  use constant_module
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8)  :: t,x1,x2 !Given time,space and velocity
  !----------------------------------------------------------------------------------
  !Calculations:
  f2=-(omega_02+2.0D0*a_0*cos(omega*t))*sin(x1)-gamma*x2
end function f2
!====================================================================================
subroutine RK(T,X1,X2,Ti,Tf,X10,X20,Nt)
  use constant_module
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  integer                :: Nt
  real(8),dimension(Nt)  :: T,X1,X2
  real(8)                :: Ti,Tf,X10,X20
  integer                :: i,j
  real(8)                :: dt         !The step of the time dt=(Tf-Ti)/(Nt-1)
  real(8)                :: TS,X1S,X2S !Local variables that hold the time,space and velocity for every step
  !----------------------------------------------------------------------------------
  !Initial values:
  dt=(Tf-Ti)/(Nt-1)
  T(1)=Ti
  X1(1)=X10
  X2(1)=X20
  TS=Ti
  X1S=X10
  X2S=X20
  !----------------------------------------------------------------------------------
  !Calculations:
  do i=2,Nt
     call RKSTEP(TS,X1S,X2S,dt) !We call a subroutine to calculate the next time, next place and next velocity
     T(i)=TS
     X1(i)=X1S
     X2(i)=X2S   !New values to our arrays
  end do
  !----------------------------------------------------------------------------------
end subroutine RK
!====================================================================================
subroutine RKSTEP(TS,X1S,X2S,dt)
  use constant_module 
  implicit none
  !----------------------------------------------------------------------------------
  !Declaration of variables:
  real(8)   :: TS,X1S,X2S,dt
  real(8)   :: f1,f2           !We declare our functions that we are going to use
  real(8)   :: k11,k12,k13,k14
  real(8)   :: k21,k22,k23,k24
  real(8)   :: h,h2,h6         !We have that h is dt, h2=h/2,h6=h/6
  !----------------------------------------------------------------------------------
  !Initial values:
  h  =  dt         !h =dt, integration step                                                                                                                                                                                               
  h2 =  0.5D0 * h  !h2=h/2                                                                                                                                                                                                                
  h6 =h/(6.0D0)    !h6=h/6
  !----------------------------------------------------------------------------------
  !Calculations:
  k11=f1(TS,X1S,X2S)
  k21=f2(TS,X1S,X2S)
  k12=f1(TS+h2,X1S+h2*k11,X2S+h2*k21)
  k22=f2(TS+h2,X1S+h2*k11,X2S+h2*k21)
  k13=f1(TS+h2,X1S+h2*k12,X2S+h2*k22)
  k23=f2(TS+h2,X1S+h2*k12,X2S+h2*k22)
  k14=f1(TS+h ,X1S+h *k13,X2S+h *k23)
  k24=f2(TS+h ,X1S+h *k13,X2S+h *k23)
  TS  =TS+h
  X1S =X1S+h6*(k11+2.0D0*(k12+k13)+k14)
  X2S =X2S+h6*(k21+2.0D0*(k22+k23)+k24) !Next steps
  !The x1 because it's an angle it must be between [-pi,pi]
  if( X1S >  PI) X1S = X1S - PI2
  if( X1S < -PI) X1S = X1S + PI2
end subroutine RKSTEP
!====================================================================================

