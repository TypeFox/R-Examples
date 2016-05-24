"logdensity.gp" <-
function(object,...){
  # calculate prior log density of the process coefficients
  m1=object$gridsize[1]
  m2=object$gridsize[2]
  # first create a vector of the real and imaginary components of the coefficients
  a1r=Re(object$coeff[1,1])
  a2r=Re(object$coeff[(m1/2+1),1])
  b1r=Re(object$coeff[2:(m1/2),1])
  b1i=Im(object$coeff[2:(m1/2),1])
  if(object$d==2){
    a3r=Re(object$coeff[1,(m2/2+1)])
    a4r=Re(object$coeff[(m1/2+1),(m2/2+1)])
    b2r=Re(object$coeff[1,2:(m2/2)])
    b3r=Re(object$coeff[(m1/2+1),2:(m2/2)])
    b4r=Re(object$coeff[2:(m1/2),(m2/2+1)])
    b2i=Im(object$coeff[1,2:(m2/2)])
    b3i=Im(object$coeff[(m1/2+1),2:(m2/2)])
    b4i=Im(object$coeff[2:(m1/2),(m2/2+1)])
    c1r=Re(object$coeff[2:(m1/2),2:(m2/2)])
    c1i=Im(object$coeff[2:(m1/2),2:(m2/2)])
    c2r=Re(object$coeff[(m1/2+2):m1,2:(m2/2)])
    c2i=Im(object$coeff[(m1/2+2):m1,2:(m2/2)])
    coeff.vec=c(a2r,a3r,a4r,b1r,b2r,b3r,b4r,b1i,b2i,b3i,b4i,c1r,c1i,c2r,c2i)
    # create vector of component variances
    variance.vec=c(object$variances[(m1/2+1),1],object$variances[1,(m2/2+1)],object$variances[(m1/2+1),(m2/2+1)],object$variances[2:(m1/2),1],object$variances[1,2:(m2/2)],object$variances[(m1/2+1),2:(m2/2)],object$variances[2:(m1/2),(m2/2+1)],object$variances[2:(m1/2),1],object$variances[1,2:(m2/2)],object$variances[(m1/2+1),2:(m2/2)],object$variances[2:(m1/2),(m2/2+1)],object$variances[2:(m1/2),2:(m2/2)],object$variances[2:(m1/2),2:(m2/2)],object$variances[(m1/2+2):m1,2:(m2/2)],object$variances[(m1/2+2):m1,2:(m2/2)])
  } else{
    coeff.vec=c(a2r,b1r,b1i)
    variance.vec=c(object$variances[(m1/2+1),1],object$variances[2:(m1/2),1],object$variances[2:(m1/2),1])
  }
  if(object$const.fixed){
    extra=0
  } else{
    extra=dnorm(a1r,0,sqrt(object$variances[1,1]),log=TRUE)
  } 
  # normal prior log density
  return(extra+sum(dnorm(coeff.vec,0,sqrt(variance.vec),log=TRUE)))
}
