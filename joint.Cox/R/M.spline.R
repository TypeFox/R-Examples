M.spline=function(time,xi1,xi3){
  if(min(time)<xi1){warning("time should be larger than xi1")}
  if(max(time)>xi3){warning("time should be smaller than xi3")}
  
  D=(xi3-xi1)/2
  xi2=xi1+D

  z1=(time-xi1)/D;z2=(time-xi2)/D;z3=(time-xi3)/D
  
  M1=-(4*z2^3/D)*(time<xi2)+0*(time>=xi2)
  M2=(7*z1^3-18*z1^2+12*z1)/2/D*(time<xi2)-z3^3/2/D*(time>=xi2)
  M3=(-2*z1^3+3*z1^2)/D*(time<xi2)+(2*z2^3-3*z2^2+1)/D*(time>=xi2)
  M4=z1^3/2/D*(time<xi2)+(-7*z2^3+3*z2^2+3*z2+1)/2/D*(time>=xi2)
  M5=4*z2^3/D*(time>=xi2)
  
  cbind(M1,M2,M3,M4,M5)
}