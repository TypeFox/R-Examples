I.spline=function(time,xi1,xi3){
  if(min(time)<xi1){warning("time should be larger than xi1")}
  if(max(time)>xi3){warning("time should be smaller than xi3")}
  
  D=(xi3-xi1)/2
  xi2=xi1+D
  
  z1=(time-xi1)/D;z2=(time-xi2)/D;z3=(time-xi3)/D

  I1=(1-z2^4)*(time<xi2)+1*(time>=xi2)
  I2=( 7/8*z1^4-3*z1^3+3*z1^2 )*(time<xi2)+( 1-z3^4/8 )*(time>=xi2)
  I3=( -z1^4/2+z1^3 )*(time<xi2)+( 1/2+z2^4/2-z2^3+z2 )*(time>=xi2)
  I4=( z1^4/8 )*(time<xi2)+( 1/8-7/8*z2^4+1/2*z2^3+3/4*z2^2+1/2*z2 )*(time>=xi2)
  I5=z2^4*(time>=xi2)
  
  cbind(I1,I2,I3,I4,I5)
}