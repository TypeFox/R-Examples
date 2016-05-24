spCor <- function(D, cov.model, phi, nu=NULL){
  if(cov.model == "exponential"){
    R <- exp(-phi*D)
  }else if(cov.model == "matern"){
    R <- (D*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=D*phi, nu=nu)
    diag(R) <- 1
  }else if(cov.model == "gaussian"){
    R <- exp(-1*((phi*D)^2))
  }else if(cov.model == "spherical"){
    R <- D
    R[TRUE] <- 1
    R[D > 0 & D < 1/phi] <- 1-1.5*phi*D[D > 0 & D <= 1/phi]+0.5*((phi*D[D > 0 & D <= 1/phi])^3)
    R[D >= 1/phi] <- 0   
  }else{
    stop("error: in spCor, specified cov.model '",cov.model,"' is not a valid option; choose, from gaussian, exponential, matern, spherical.")
  }
  R
}
