#'  Circular variance and concentration parameter 
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the Cayley distribution.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @@cite mardia2000
#'  @seealso \code{\link{Cayley}}
#'  @export
#'  @examples
#'  #Find the concentration parameter for circular variances 0.25, 0.5, 0.75
#'  cayley.kappa(0.25)
#'  cayley.kappa(0.5)
#'  cayley.kappa(0.75)

cayley.kappa<-function(nu){
  (3/nu)-2
}

fisher.nu.kappa<-function(kappa,nu){
  ((1-(besselI(2*kappa,1)-.5*besselI(2*kappa,2)-.5*besselI(2*kappa,0))/(besselI(2*kappa,0)-besselI(2*kappa,1)))-nu)^2
}

#'  Circular variance and concentration parameter 
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the matrix-Fisher distribution.  For numerical stability, a maximum \eqn{\kappa} of 350 is returned.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @seealso \code{\link{Fisher}}
#'  @@cite mardia2000
#'  @export
#'  @examples
#'  #Find the concentration parameter for circular variances 0.25, 0.5, 0.75
#'  fisher.kappa(0.25)
#'  fisher.kappa(0.5)
#'  fisher.kappa(0.75)

fisher.kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(fisher.nu.kappa,interval=c(0,350),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}


mises.nu.kappa<-function(kappa,nu){
  (1-besselI(kappa,1)/besselI(kappa,0)-nu)^2
}

#'  Circular variance and concentration parameter
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the circular-von Mises distribution.  For numerical stability, a maximum \eqn{\kappa} of 500 is returned.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @seealso \code{\link{Mises}}
#'  @@cite mardia2000
#'  @export
#'  @examples
#'  #Find the concentration parameter for circular variances 0.25, 0.5, 0.75
#'  vmises.kappa(0.25)
#'  vmises.kappa(0.5)
#'  vmises.kappa(0.75)

vmises.kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(mises.nu.kappa,interval=c(0,500),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}


maxwell.nu.kappa <- function(kappa,nu){
  ((1-(2*sqrt(kappa*pi)*exp(-kappa*pi*pi)+(2*kappa-1)*exp(-1/(4*kappa))/(2*kappa)))-nu)^2
}


#'  Circular variance and concentration parameter
#'  
#'  Return the concentration parameter that corresponds to a given circular variance.
#'  
#'  The concentration parameter \eqn{\kappa} does not translate across circular distributions.  A commonly used
#'  measure of spread in circular distributions that does translate
#'  is the circular variance defined as \eqn{\nu=1-E[\cos(r)]}{\nu=1-E[cos(r)]} where \eqn{E[\cos(r)]}{E[cos(r)]} 
#'  is the mean resultant length.  See \cite{mardia2000} for more details.
#'  This function translates the circular variance \eqn{\nu} into the corresponding concentration parameter \eqn{\kappa}
#'  for the modified Maxwell-Boltzmann distribution.  For numerical stability, a maximum \eqn{\kappa} of 1000 is returned.
#'  
#'  @param nu circular variance
#'  @return Concentration parameter corresponding to nu.
#'  @seealso \code{\link{Maxwell}}
#'  @export
#'  @examples
#'  #Find the concentration parameter for circular variances 0.25, 0.5, 0.75
#'  maxwell.kappa(0.25)
#'  maxwell.kappa(0.5)
#'  maxwell.kappa(0.75)

maxwell.kappa<-function(nu){
  
  kappa<-rep(0,length(nu))
  
  for(i in 1:length(nu))
    kappa[i]<-optimize(maxwell.nu.kappa,interval=c(0,1000),tol=.00001,nu=nu[i])$minimum
  
  return(kappa)
}
