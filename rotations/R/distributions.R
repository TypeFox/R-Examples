dmkern <- function(r,kappa){
  return(2*kappa*sqrt(kappa/pi)*(r^2)*exp(-kappa*(r^2)))
}

#' Angular distributions
#' 
#' Density, distribution function and random variate generation for symmetric probability distributions in the rotations package.
#' 
#' The functions for the density function and random variate generation are named in the usual form dxxxx, pxxxx and rxxxx, 
#' respectively.        
#' \itemize{
#' 	\item See \code{\link{Cayley}} for the Cayley distribution.
#' 	\item See \code{\link{Fisher}} for the matrix Fisher distribution.
#' 	\item See \code{\link{Haar}} for the uniform distribution on the circle.
#' 	\item See \code{\link{Maxwell}} for the Maxwell-Boltzmann distribution on the circle.
#' 	\item See \code{\link{Mises}} for the von Mises-Fisher distribution.
#' }
#' 
#' @name Angular-distributions

NULL



arsample <- function(f, g, M, kappa, Haar, ...) {
  #generate a random observation from target density f
  found = FALSE
  while (!found) {
    x <- g(1, ...)
    y <- runif(1, min = 0, max = M)
    if (y < f(x, kappa, Haar)) 
      found = TRUE
  }
  return(x)
  # arsample(f, g, M, kappa, ...)
}



arsample.unif <- function(f, M, ...) {
  #generate a random observation from target density f assuming g is uniform
  found = FALSE
  while (!found) {
    x <- runif(1, -pi, pi)
    y <- runif(1, min = 0, max = M)
    if (y < f(x, ...)) 
      found = TRUE
  }
  return(x)
  # arsample.unif(f, M, ...)
}


rar <- function(n, f, M, ...) {
  res <- vector("numeric", length = n)
  for (i in 1:n) res[i] <- arsample.unif(f, M, ...)
  return(res)
}


#' The symmetric Cayley distribution
#'
#' Density, distribution function and random generation for the Cayley distribution with concentration \code{kappa} \eqn{\kappa}.
#'
#' The symmetric Cayley distribution with concentration \eqn{\kappa} has density 
#' \deqn{C_C(r |\kappa)=\frac{1}{\sqrt{\pi}} \frac{\Gamma(\kappa+2)}{\Gamma(\kappa+1/2)}2^{-(\kappa+1)}(1+\cos r)^\kappa(1-\cos r).}{C(r |\kappa)= \Gamma(\kappa+2)(1+cos r)^\kappa(1-cos r)/[\Gamma(\kappa+1/2)2^(\kappa+1)\sqrt\pi].}
#' The Cayley distribution is equivalent to the de la Vallee Poussin distribution of \cite{Schaeben1997}.
#'
#' @name Cayley
#' @aliases Cayley rcayley dcayley
#' @param r,q vector of quantiles.
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required.
#' @param kappa concentration parameter.
#' @param nu circular variance, can be used in place of \code{kappa}.
#' @param Haar logical; if TRUE density is evaluated with respect to the Haar measure.
#' @param lower.tail logical; if TRUE (default) probabilities are \eqn{P(X\leq x)}{P(X\le x)} otherwise, \eqn{P(X>x)}.
#' @return  \item{dcayley}{gives the density}
#'          \item{pcayley}{gives the distribution function}
#'          \item{rcayley}{generates a vector of random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package.
#' @@cite Schaeben1997 leon2006
#' @examples
#' r <- seq(-pi, pi, length = 500)
#' 
#' #Visualize the Cayley density fucntion with respect to the Haar measure
#' plot(r, dcayley(r, kappa = 10), type = "l", ylab = "f(r)")
#' 
#' #Visualize the Cayley density fucntion with respect to the Lebesgue measure
#' plot(r, dcayley(r, kappa = 10, Haar = FALSE), type = "l", ylab = "f(r)")
#' 
#' #Plot the Cayley CDF
#' plot(r,pcayley(r,kappa = 10), type = "l", ylab = "F(r)")
#' 
#' #Generate random observations from Cayley distribution
#' rs <- rcayley(20, kappa = 1)
#' hist(rs, breaks = 10)

NULL


#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

dcayley <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- cayley.kappa(nu)
  
 	den <- 0.5 * gamma(kappa + 2)/(sqrt(pi) * 2^kappa * gamma(kappa + 0.5)) * (1 + cos(r))^kappa * (1 - cos(r))
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
}

#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

pcayley<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  if(!is.null(nu))
    kappa <- cayley.kappa(nu)
  
  n<-length(q)
  cdf<-rep(NA,n)
  a<-(2*kappa+1)*(1-cos(q))/(3*(1+cos(q)))
  pos<-which(q>0)
  
  cdf<-pf(a,3,2*kappa+1)
  cdf[pos]<-cdf[pos]/2+0.5
  cdf[-pos]<-0.5-cdf[-pos]/2

  if(lower.tail) 
    return(cdf) else return((1-cdf))
}


#' @rdname Cayley
#' @aliases Cayley dcayley pcayley rcayley
#' @export

rcayley <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- cayley.kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  #This has also been coded in C++ but the gains aren't worth the pains so
  #I will keep it in R for now
  bet <- rbeta(n, kappa + 0.5, 3/2)
  theta <- acos(2 * bet - 1) * (1 - 2 * rbinom(n, 1, 0.5))
  return(theta)
}

#' The matrix-Fisher distribution
#'
#' Density, distribution function and random generation for the matrix-Fisher distribution with concentration \code{kappa} \eqn{\kappa}.
#'
#' The matrix-Fisher distribution with concentration \eqn{\kappa} has density
#' \deqn{C_\mathrm{{F}}(r|\kappa)=\frac{1}{2\pi[\mathrm{I_0}(2\kappa)-\mathrm{I_1}(2\kappa)]}e^{2\kappa\cos(r)}[1-\cos(r)]}{C(r|\kappa)=exp[2\kappa cos(r)][1-cos(r)]/(2\pi[I0(2\kappa)-I1(2\kappa)])}
#' with respect to Lebesgue measure where \eqn{\mathrm{I}_p(\cdot)}{Ip()} denotes the Bessel function of order \eqn{p} defined as  
#' \eqn{\mathrm{I}_p(\kappa)=\frac{1}{2\pi}\int_{-\pi}^{\pi}\cos(pr)e^{\kappa\cos r}dr}{Ip(\kappa)}.  If \code{kappa>354} then random deviates
#' are generated from the \code{\link{Cayley}} distribution because they agree closely for large \code{kappa} and generation is 
#' more stable from the Cayley distribution.
#' 
#' For large \eqn{\kappa}, the Bessel functon gives errors so a large \eqn{\kappa} approximation to the matrix-Fisher
#' distribution is used instead, which is the Maxwell-Boltzmann density.
#'
#' @name Fisher
#' @aliases Fisher dfisher rfisher pfisher
#' @param r,q vector of quantiles.
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required.
#' @param kappa concentration parameter.
#' @param nu circular variance, can be used in place of \code{kappa}.
#' @param Haar logical; if TRUE density is evaluated with respect to the Haar measure.
#' @param lower.tail  logical; if TRUE (default), probabilities are \eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.
#' @return  \item{dfisher}{gives the density}
#'          \item{pfisher}{gives the distribution function}
#'          \item{rfisher}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package.
#' @examples
#' r <- seq(-pi, pi, length = 500)
#' 
#' #Visualize the matrix Fisher density fucntion with respect to the Haar measure
#' plot(r, dfisher(r, kappa = 10), type = "l", ylab = "f(r)")
#' 
#' #Visualize the matrix Fisher density fucntion with respect to the Lebesgue measure
#' plot(r, dfisher(r, kappa = 10, Haar = FALSE), type = "l", ylab = "f(r)")
#' 
#' #Plot the matrix Fisher CDF
#' plot(r,pfisher(r,kappa = 10), type = "l", ylab = "F(r)")
#' 
#' #Generate random observations from matrix Fisher distribution
#' rs <- rfisher(20, kappa = 1)
#' hist(rs, breaks = 10)

NULL

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export

dfisher <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- fisher.kappa(nu)
  
  n<-length(r)
  den<-rep(0,n)
  if(kappa<200){
    #For small kappa use the matrix Fisher density directly
 	  den <- exp(2 * kappa * cos(r)) * (1 - cos(r))/(2 * pi * (besselI(2 * kappa, 0) - besselI(2 * kappa, 1)))
  }
  else{
    #besselI function exhibits bad behavior for large kappa so use the Maxwell-Boltzman approx then
    den <- 2*kappa*sqrt(kappa/pi)*r^2*exp(-kappa*r^2)
  }
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
  
}

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export

pfisher<-function(q,kappa=1, nu=NULL, lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dfisher,-pi,q[i],kappa,nu,Haar=FALSE)$value,1),0)
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @rdname Fisher
#' @aliases Fisher dfisher pfisher rfisher
#' @export


rfisher <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises.kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
    n<-lenn
  
  theta<-rfisherCpp(n,kappa)
  return(theta)
}

#' Uniform distribution
#'
#' Density, distribution function and random generation for the uniform distribution.
#'
#' The uniform distribution
#' has density \deqn{C_U(r)=\frac{[1-cos(r)]}{2\pi}}{C(r)=[1-cos(r)]/2\pi} with respect to the Lebesgue
#' measure.  The Haar measure is the volume invariant measure for SO(3) that plays the role
#' of the uniform measure on SO(3) and \eqn{C_U(r)}{C(r)} is the angular distribution that corresponds
#' to the uniform distribution on SO(3), see \code{\link{UARS}}.  The uniform distribution with respect to the Haar measure is given
#' by \deqn{C_U(r)=\frac{1}{2\pi}.}{C(r)=1/(2\pi).}  Because the uniform distribution
#' with respect to the Haar measure gives a horizontal line at 1 
#' with respect to the Lebesgue measure, we called this distribution 'Haar.'
#'
#' @name Haar
#' @aliases Haar dhaar phaar rhaar
#' @param r,q vector of quantiles.
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required.
#' @param lower.tail  logical; if TRUE (default), probabilities are \eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.
#' @return  \item{dhaar}{gives the density}
#'          \item{phaar}{gives the distribution function}
#'          \item{rhaar}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package.
#' @examples
#' r <- seq(-pi, pi, length = 1000)
#' 
#' #Visualize the uniform distribution with respect to Lebesgue measure
#' plot(r, dhaar(r), type = "l", ylab = "f(r)")
#' 
#' #Visualize the uniform distribution with respect to Haar measure, which is
#' #a horizontal line at 1
#' plot(r, 2*pi*dhaar(r)/(1-cos(r)), type = "l", ylab = "f(r)")
#' 
#' #Plot the uniform CDF
#' plot(r,phaar(r), type = "l", ylab = "F(r)")
#' 
#' #Generate random observations from uniform distribution
#' rs <- rhaar(50)
#' 
#' #Visualize on the real line
#' hist(rs, breaks = 10)
#' 


NULL

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

dhaar <- function(r){
	
	den <-(1 - cos(r))/(2 * pi)
	
	return(den)
} 

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

phaar<-function(q,lower.tail=TRUE){
  
  cdf<-(q-sin(q)+pi)/(2*pi)
  
  ind<-which(cdf>1)
  cdf[ind]<-1
  
  ind2<-which(cdf<0)
  cdf[ind2]<-0
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @rdname Haar
#' @aliases Haar dhaar phaar rhaar
#' @export

rhaar<-function(n){
	
	lenn<-length(n)
	if(lenn>1)
		n<-lenn
	
  return(rar(n, dhaar, 1/pi))
}

#' The modified Maxwell-Boltzmann distribution
#'
#' Density, distribution function and random generation for the Maxwell-Boltzmann distribution with 
#' concentration \code{kappa} \eqn{\kappa} restricted to the range \eqn{[-\pi,\pi)}.
#'
#' The Maxwell-Boltzmann distribution with concentration \eqn{\kappa} has density
#' \deqn{C_\mathrm{{M}}(r|\kappa)=2\kappa\sqrt{\frac{\kappa}{\pi}}r^2e^{-\kappa r^2}}{C(r|\kappa)=2\kappa(\kappa/\pi)^(1/2)r^2exp(-\kappa r^2)}
#' with respect to Lebesgue measure.  The usual expression for the Maxwell-Boltzmann distribution can be recovered by
#' setting \eqn{a=(2\kappa)^0.5}.
#'
#' @name Maxwell
#' @aliases Maxwell rmaxwell dmaxwell pmaxwell
#' @param r,q vector of quantiles.
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required.
#' @param kappa concentration parameter.
#' @param nu circular variance, can be used in place of \code{kappa}.
#' @param Haar logical; if TRUE density is evaluated with respect to the Haar measure.
#' @param lower.tail logical; if TRUE (default) probabilities are \eqn{P(X\leq x)}{P(X\le x)} otherwise, \eqn{P(X>x)}.
#' @return  \item{dmaxwell}{gives the density}
#'          \item{pmaxwell}{gives the distribution function}
#'          \item{rmaxwell}{generates a vector of random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package.
#' @@cite bingham2010
#' @examples
#' r <- seq(-pi, pi, length = 500)
#' 
#' #Visualize the Maxwell-Boltzmann density fucntion with respect to the Haar measure
#' plot(r, dmaxwell(r, kappa = 10), type = "l", ylab = "f(r)")
#' 
#' #Visualize the Maxwell-Boltzmann density fucntion with respect to the Lebesgue measure
#' plot(r, dmaxwell(r, kappa = 10, Haar = FALSE), type = "l", ylab = "f(r)")
#' 
#' #Plot the Maxwell-Boltzmann CDF
#' plot(r,pmaxwell(r,kappa = 10), type = "l", ylab = "F(r)")
#' 
#' #Generate random observations from Maxwell-Boltzmann distribution
#' rs <- rmaxwell(20, kappa = 1)
#' hist(rs, breaks = 10)

NULL


#' @name Maxwell
#' @aliases Maxwell rmaxwell dmaxwell pmaxwell
#' @export

dmaxwell <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- maxwell.kappa(nu)
  
  den <- dmkern(r = r, kappa = kappa)
  
  if(kappa<1.9){
    den <- den/(1-2*integrate(dmkern,lower = pi, upper = Inf,kappa = kappa)$value)
  }
  
  zeros <- which(abs(r)>pi)
  if(length(zeros)>0)
    den[zeros] <- 0
  
  if (Haar) 
    return(den/(1 - cos(r))) else return(den)
}


#' @name Maxwell
#' @aliases Maxwell rmaxwell dmaxwell pmaxwell
#' @export

pmaxwell<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dmaxwell,-pi,q[i],kappa,nu,Haar=FALSE)$value,1),0)
  
  if(lower.tail)
    return(cdf) else return((1-cdf))
}

#' @name Maxwell
#' @aliases Maxwell rmaxwell dmaxwell pmaxwell
#' @export

rmaxwell <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- maxwell.kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
    n<-lenn
  
  theta<-rmbCpp(n,kappa)
  return(theta)
}

#' The circular-von Mises distribution
#'
#' Density, distribution function and random generation for the circular-von Mises distribution with concentration \code{kappa} \eqn{\kappa}.
#' 
#' The circular von Mises distribution with concentration \eqn{\kappa} has density
#' \deqn{C_\mathrm{M}(r|\kappa)=\frac{1}{2\pi \mathrm{I_0}(\kappa)}e^{\kappa cos(r)}.}{C(r|\kappa)=exp[\kappa cos(r)]/[2\pi I(\kappa)]}
#' where \eqn{\mathrm{I_0}(\kappa)}{I(\kappa)} is the modified Bessel function of order 0.
#'
#' @name Mises
#' @aliases Mises dvmises rvmises pvmises
#' @param r,q vector of quantiles
#' @param n number of observations.  If \code{length(n)>1}, the length is taken to be the number required.
#' @param kappa concentration parameter.
#' @param nu circular variance, can be used in place of \code{kappa}.
#' @param Haar logical; if TRUE density is evaluated with respect to the Haar measure.
#' @param lower.tail  logical; if TRUE (default), probabilities are \eqn{P(X \le x)} otherwise, \eqn{P(X > x)}.
#' @return  \item{dvmises}{gives the density}
#'          \item{pvmises}{gives the distribution function}
#'          \item{rvmises}{generates random deviates}
#' @seealso \link{Angular-distributions} for other distributions in the rotations package.
#' @examples
#' r <- seq(-pi, pi, length = 500)
#' 
#' #Visualize the von Mises density fucntion with respect to the Haar measure
#' plot(r, dvmises(r, kappa = 10), type = "l", ylab = "f(r)", ylim = c(0, 100))
#' 
#' #Visualize the von Mises density fucntion with respect to the Lebesgue measure
#' plot(r, dvmises(r, kappa = 10, Haar = FALSE), type = "l", ylab = "f(r)")
#'
#' #Plot the von Mises CDF
#' plot(r,pvmises(r,kappa = 10), type = "l", ylab = "F(r)")
#'   
#' #Generate random observations from von Mises distribution
#' rs <- rvmises(20, kappa = 1)
#' hist(rs, breaks = 10)

NULL

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

dvmises <- function(r, kappa = 1, nu = NULL, Haar = TRUE) {
  
  if(!is.null(nu))
    kappa <- vmises.kappa(nu)
  
  den <- 1/(2 * pi * besselI(kappa, 0)) * exp(kappa * cos(r))
  
  #if(!lower.tail)
  #	den<-1-den
  
  if (Haar) {
    return(den/(1 - cos(r)))
  } else {
    return(den)
  }
}

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

pvmises<-function(q,kappa=1,nu=NULL,lower.tail=TRUE){
  
  n<-length(q)
  cdf<-rep(NA,n)
  
  for(i in 1:n)
    cdf[i]<-max(min(integrate(dvmises,-pi,q[i],kappa,nu,Haar=FALSE)$value,1),0)
  
  if(lower.tail) 
    return(cdf) else return((1-cdf))
}

#' @rdname Mises
#' @aliases Mises dvmises pvmises rvmises
#' @export

rvmises <- function(n, kappa = 1, nu = NULL) {
  
  if(!is.null(nu))
    kappa <- vmises.kappa(nu)
  
  lenn<-length(n)
  if(lenn>1)
  	n<-lenn
  
  theta<-rvmisesCPP(n,kappa)
  return(theta)
}



#' Generic UARS Distribution
#'
#' Density, distribution function and random generation for the the generic uniform axis-random spin (UARS) class of distributions.
#' 
#' For the rotation R with central orientation S and concentration \eqn{\kappa} the UARS density is given by 
#' \deqn{f(R|S,\kappa)=\frac{4\pi}{3-tr(S^\top R)}C(\cos^{-1}[tr(S^\top R)-1]/2|\kappa)}{f(R|S,\kappa)=4\pi C(acos[tr(S'R)-1]/2|\kappa)/[3-tr(S'R)]}
#' where \eqn{C(r|\kappa)} is one of the \link{Angular-distributions}.
#'
#' @name UARS
#' @aliases UARS puars duars ruars
#' @param R Value at which to evaluate the UARS density.
#' @param n number of observations. If \code{length(n)>1}, the length is taken to be the number required.
#' @param dangle The function to evaluate the angles from, e.g. dcayley, dvmises, dfisher, dhaar.
#' @param pangle The form of the angular density, e.g. pcayley, pvmises, pfisher, phaar.
#' @param rangle The function from which to simulate angles, e.g. rcayley, rvmises, rhaar, rfisher.
#' @param S central orientation of the distribution.
#' @param kappa concentration parameter.
#' @param space indicates the desired representation: matrix ("SO3") or quaternion ("Q4").
#' @param ... additional arguments.
#' @return  \item{duars}{gives the density}
#'          \item{puars}{gives the distribution function.  If pangle is left empty, the empirical CDF is returned.}
#'          \item{ruars}{generates random deviates}
#' @seealso For more on the angular distribution options see \link{Angular-distributions}.
#' @@cite bingham09
#' @examples
#' #Generate random rotations from the Cayley-UARS distribution with central orientation 
#' #rotated about the y-axis through pi/2 radians
#' S <- as.SO3(c(0, 1, 0), pi/2)
#' Rs <- ruars(20, rangle = rcayley, kappa = 1, S = S)
#' 
#' rs <- mis.angle(Rs-S)                          #Find the associated misorientation angles
#' frs <- duars(Rs, dcayley, kappa = 10, S = S)   #Compute UARS density evaluated at each rotations
#' plot(rs, frs) 
#' 
#' cdf <- puars(Rs, pcayley, S = S)               #By supplying 'pcayley', it is used to compute the
#' plot(rs, cdf)                                  #the CDF
#' 
#' ecdf <- puars(Rs, S = S)                       #No 'puars' arguement is supplied so the empirical
#' plot(rs, ecdf)                                 #cdf is returned

NULL


#' @rdname UARS
#' @aliases UARS duars puars ruars
#' @export

duars<-function(R,dangle,S=id.SO3,kappa=1,...){
	
	R<-formatSO3(R)
	rs<-mis.angle(R-S)
	cr<-dangle(rs,kappa,...)	
	trStO<-2*cos(rs)+1
	
	den<-4*pi*cr/(3-trStO)
	
	return(den)
}

#' @rdname UARS
#' @aliases UARS duars puars ruars
#' @export

puars<-function(R,pangle=NULL,S=id.SO3,kappa=1,...){
	
	#This is not a true CDF, but it will work for now
	R<-formatSO3(R)
	rs<-mis.angle(R-S)
  
	if(is.null(pangle)){
		
		n<-length(rs)
		cr<-rep(0,n)
		
		for(i in 1:length(rs))
			cr[i]<-length(which(rs<=rs[i]))/n
		
	}else{		
		cr<-2*(pangle(rs,kappa,...)-.5)
	}
	
  
	
	return(cr)
	
}

#' @rdname UARS
#' @aliases UARS duars puars ruars
#' @export

ruars<-function(n,rangle,S=NULL,kappa=1,space="SO3",...){
  
  r<-rangle(n,kappa,...)
  Rs<-genR(r,S,space)
  
  return(Rs)
}
