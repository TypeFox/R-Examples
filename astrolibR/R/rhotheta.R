rhotheta = function(p,t,e,a,i,omega,omega2,t2) {

                                        # see chapter 55
  n=360.0/p
  m=n*(t2-t)
  m=m/360.0*2.0*pi # convert m to radians

                                        # solution of kepler equation, see chapter 29, 3rd method
  f= ifelse(m > 0, 1, -1)
  m=abs(m)/2.0/pi
  m=(m-floor(m))*2.0*pi*f
  if (m < 0.0)  m=m+2.0*pi
  f=1.0
  if (m > pi) f=-1.0
  if (m > pi) m=2.0*pi-m
  e0=pi/2.0
  d=pi/4.0
  for (j in 1:33) {
    m1=e0-e*sin(e0)
    sgn_m = ifelse(m-m1 > 0, 1, -1)
    e0=e0+d*sgn_m
    d=d/2.0
  }
  e0=e0*f

                                        # return to chapter 55
  r=a*(1.0-e*cos(e0))
  nu=2.0*atan(sqrt((1.0+e)/(1.0-e))*tan(e0/2.0))
  my_omega2=omega2/180.0*pi # convert variables in radians and copy them to a new variable to prevent changes to the input parameter
  my_i=i/180.0*pi
  my_omega=omega/180.0*pi
  theta=my_omega+atan2(sin(nu+my_omega2)*cos(my_i),cos(nu+my_omega2))
  rho=r*cos(nu+my_omega2)/cos(theta-my_omega)
  theta=theta*180.0/pi # convert theta to degree

  theta = theta %% 360 # force theta to be in 0..360 range
  list(rho=rho, theta=theta)
}
