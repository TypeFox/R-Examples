helio_rv = function(hjd,t,p,v0,k,e,omega,maxiter=100) {

  if(missing(omega) && missing(e) ){
    e = 0.0
    omega = 0.0
  }
  dtor = pi/180
  m=2*pi*( ((hjd-t)/p) %% 1)

  e1=m + e*sin(m)  + ((e^2)*sin(2.0*m)/2.0)
  
  for(i in 1:maxiter) {
    e0=e1 
    m0 = e0 - e*sin(e0)
    e1 = e0 + (m-m0)/(1.0 - e*cos(e0))
    tmp = e1!=0
    if(sum(tmp)>0) {
      test = max(abs( (e1[tmp]-e0[tmp])/e1[tmp]))
      if(test<1e-8) break
    }
  }
  if(i==maxiter) return(list(status=paste("divergence after",maxiter,"iterations"),rv=0))
  nu=2*atan(sqrt((1.e0 + e)/(1 - e))*tan(e1/2))
  rv = (k*(cos(nu+dtor*omega) + (e*cos(dtor*omega))))+v0
  return(list(status=paste("convergence after",i,"iterations"),rv=rv))
}
