##' semi-infinite DBR
##'
##' periodic structure of dielectric layers
##' @title dbr_analytic
##' @export
##' @param wavelength in nm
##' @param lambda0 central wavelength of the stopband
##' @param n1 real refractive index for odd layers
##' @param n2 real refractive index for even layers
##' @param d1 odd layer thickness in nm
##' @param d2 even layer thickness in nm
##' @param nleft real refractive index for incident medium
##' @param ... ignored
##' @return data.frame with complex reflectivity
##' @references Amir and Vukusic, 2013, arXiv:1209.3776v2
##' @note issue at lambda0/2 needs investigating
##' @author baptiste Auguie
##' @family dbr user_level
dbr_analytic <- function(wavelength, lambda0,
                         n1, n2, nleft, 
                         d1=lambda0/4/n1,
                         d2=lambda0/4/n2, ...){
  k0 <- 2*pi/wavelength
  rlr <- (nleft - n1)/(nleft + n1)
  rrl <- exp(1i*k0*n1*d1)*(n1 - nleft)/(nleft + n1)
  tlr <- 2*nleft*exp(1i*k0*n1*d1/2)/(nleft + n1)
  trl <- 2*n1*exp(1i*k0*n1*d1/2)/(nleft + n1)
  r12 <- (n1-n2)/(n1+n2)
  t12 <- 2*n1/(n1+n2)
  t21 <- 2*n2/(n1+n2)
  r21 <- -r12
  r <- r12*exp(1i*k0*n1*d1) +(r21*t12*t21*exp(1i*k0*n1*d1+2i*k0*n2*d2)) / 
    (1 - (r21*exp(1i*k0*n2*d2))^2)
  rfn <- -(r+Conj(r)) / (2*Mod(r)^2) - sqrt((r+Conj(r))^2/(4*Mod(r)^4)-1)
  rfp <- -(r+Conj(r)) / (2*Mod(r)^2) + sqrt((r+Conj(r))^2/(4*Mod(r)^4)-1)
  rtot1 <- rlr + (rfn*tlr*trl) / (1 - rrl*rfn)
  rtot2 <- rlr + (rfp*tlr*trl) / (1 - rrl*rfp)
  
  rtot <- ifelse(wavelength<lambda0, rtot1, -Conj(rtot2))
  
  data.frame(wavelength=wavelength, r=rtot, R=Mod(rtot)^2)
}

##' combine layer
##'
##' reflection coefficient for a layer
##' @title combine_layer
##' @export
##' @param r1 reflection coefficient
##' @param r2 reflection coefficient
##' @param kd k*d
##' @return combined complex reflectivity
##' @author baptiste Auguie
##' @family user_level
combine_layer <- function(r1, r2, kd){
  (r1+r2*exp(2i*kd)) / (1+r1*r2*exp(2i*kd))
}
