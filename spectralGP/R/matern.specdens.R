"matern.specdens" <-
function(omega,param,d=2){
  # spectral density for Matern correlation function
  # param should be in the form of (correlation scale, differentiability) (often called rho,nu)
  # this uses the form of the Matern in which rho and nu interact minimally, namely,
  # R(\tau)=\frac{1}{\Gamma(\nu) 2^(\nu-1)} (\frac{2\sqrt{nu} tau}{\rho})^\nu K_{\nu}(\frac{2\sqrt{nu} tau}{\rho}) - see Stein (1999; p. 50), Interpolation of Spatial Data: Some Theory for Kriging
  if(is.vector(omega) || (is.matrix(omega) && min(dim(omega))==1) || (is.array(omega) && length(dim(omega))==1)){
    if(d==2 && length(omega)!=2){
      warning("You specified d=2 but gave a vector of frequencies; interpreting as d=1.\n")
      d=1
      omega=matrix(c(omega),ncol=1)
    }
    if(d==2 && length(omega)==2){
      warning("Assuming a single two-dimensional frequency.\n")
      omega=matrix(omega,nrow=1,ncol=d)
    }
  } else{
    if(!(
         (is.matrix(omega) && ncol(omega)==2) ||
         (is.array(omega) && length(dim(omega))==2 && dim(omega)[2]==2) ||
         (is.data.frame(omega) && ncol(omega)==2)
         )){
      stop("omega must be a vector or a two-dimensional matrix/array")
    }
  }
  if(is.matrix(omega) && min(dim(omega))>1 && d==1){
    stop("omega is a matrix with more than one row, but you specified d=1")
  }
  if(length(param)!=2 || param[1]<=0 || param[2]<=0){
    stop("Matern spectral density requires two positive parameters")
  }
  rho=param[1]
  nu=param[2]
  if(d==1){
    omegaTomega=omega*omega
  } else{
    omegaTomega=omega[,1]*omega[,1]+omega[,2]*omega[,2]
  }
  exp(lgamma(nu+d/2)+nu*log(4*nu)-(d/2)*log(pi)-lgamma(nu)-2*nu*log(rho*pi)-(nu+d/2)*log(4*nu/(pi*pi*rho*rho)+omegaTomega))
}
