lsf_rotate=function( deltav, vsini, epsilon=0.6){

  e1 = 2.0*(1.0 - epsilon)
  e2 = pi*epsilon/2.0
  e3 = pi*(1.0 - epsilon/3.0)
  

  npts = ceiling(2*vsini/deltav)

  if(npts %% 2==0 )npts = npts +1
  nwid = npts%/%2
  cat('nwid=',nwid)
  x = (0:(npts-1)- nwid)
  x = x*deltav/vsini  

  x1 = abs(1.0 - x^2)
  return((e1*sqrt(x1) + e2*x1)/e3)
  
}   
