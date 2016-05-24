#' @name zambia_graph
#' @title Prior precision  matrix for spatial variable in Zambia data set
#' @description This is a 57x57 matrix containing row- and columwise the regions of Zambia, and the entries
#'  define the neighbourhoodstructure. The corresponding map sambia.bnd can be downloaded from \url{http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/daten_e.html}.
#'  from the bnd file the prior precision matrix is obtained by
#'  library(BayesX)
#'  map <- read.bnd("zambia.bnd")
#'  K <- bnd2gra(map)
#' @docType data
NULL

#' @name zambia_height92
#' @title Malnutrition in Zambia
#' @description The primary goal of a statistical analysis is to determine the effect of certain socioeconomic
#'              variables of the child, the mother, and the household on the child`s nutritional
#'              condition
#'              \itemize{
#'                \item zscore child`s Z-score
#'                \item c _breastf duration of breastfeeding in months
#'                \item c_age child`s age in months
#'                \item m_agebirth mother`s age at birth in years
#'                \item m_height mother`s height in centimeter
#'                \item m_bmi mother`s body mass index
#'                \item m_education mother`s level of education
#'                \item m_work mother`s work status
#'                \item region region of residence in Zambia
#'                \item district district of residence in Zambia
#'              }
#' @format A data frame with 4421 rows and 21 variables
#' @source \url{http://www.stat.uni-muenchen.de/~kneib/regressionsbuch/daten_e.html}
#' @docType data

NULL

#' Compute Density Function of Approximated (Differentiably) Uniform Distribution. 

#' 
#' @param x  denotes the argument of the density function.  
#' @param scale the scale parameter originally defining the upper bound of the uniform distribution.
#' @param tildec denotes the ratio between scale parameter \eqn{\theta} and \eqn{s}. The latter is responsible for the closeness of
#'         the approximation to the uniform distribution. See also below for further details and the default value. 
#' @return the density.
#' @author Nadja Klein
#' @details The density of the uniform distribution for \eqn{\tau} is approximated by
#'          \deqn{p(\tau)=(1/(1+exp(\tau\tilde{c}/\theta-\tilde{c})))/(\theta(1+log(1+exp(-\tilde{c}))))}.
#'			This results in \deqn{p(\tau^2)=0.5*(\tau^2)^(-1/2)(1/(1+exp((\tau^2)^(1/2)\tilde{c}/\theta-\tilde{c})))/(\theta(1+log(1+exp(-\tilde{c}))))} for \eqn{tau^2}.
#' 			\eqn{\tilde{c}} is chosen such that \eqn{P(\tau<=\theta)>=0.95}.
#'
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#'
#' @seealso \code{\link{rapprox_unif}},\code{\link{papprox_unif}}
#' @export

dapprox_unif <- function(x,scale,tildec=13.86294) 
  {
  arg <- x^(0.5)*tildec/scale - tildec
  ret <- 0.5*x^(-0.5)*(1/(1+exp(arg)))/(scale*(1+log(1+exp(-tildec))/tildec))
  return(ret)
  }
    

NULL

#' Draw Random Numbers from Approximated (Differentiably) Uniform Distribution. 

#' 
#' @param n  number of draws.  
#' @param scale the scale parameter originally defining the upper bound of the uniform distribution.
#' @param tildec denotes the ratio between scale parameter \eqn{\theta} and \eqn{s}. The latter is responsible for the closeness of
#'         the approximation to the uniform distribution. See also below for further details and the default value. 
#' @param seed denotes the seed
#' @return n draws with density \code{\link{papprox_unif}}.
#' @author Nadja Klein
#' @details The method is based on the inversion method and the quantile function is computed numerically using \code{\link{uniroot}}.
#'
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#' @seealso \code{\link{rapprox_unif}},\code{\link{papprox_unif}}
#' @export

rapprox_unif<-function (n=100, scale,tildec=13.86294,seed=123) 
  {
  set.seed(seed)
  invF <- function(x,u)
    {
	arg <- x^(0.5)*tildec/scale 
    ( tildec*x^(0.5) -  scale*log(exp(arg)+exp(tildec))) / ( scale * (tildec + log(1+exp(-tildec)) ))+1-u
	}
  us <- runif(n)
    vec <- rep(NA,length=n)
  for(i in 1:length(us))
    {
	vec[i] <- uniroot(invF,lower=0,upper=10*scale^2,u=us[i])$root
	}
  return(vec)	
}


NULL

#' Compute Cumulative Distribution Function of Approximated (Differentiably) Uniform Distribution. 

#' 
#' @param x  denotes the argument of cumulative distribution function  
#' @param scale the scale parameter originally defining the upper bound of the uniform distribution.
#' @param tildec denotes the ratio between scale parameter \eqn{\theta} and \eqn{s}. The latter is responsible for the closeness of
#'         the approximation to the uniform distribution. See also below for further details and the default value. 
#' @return the cumulative distribution function.
#' @author Nadja Klein
#' @details The cumulative distribution function of \code{\link{dapprox_unif}} is given by
#'          \deqn{(1/(log(1+exp(-\tilde{c}))+\tilde{c}))*(\tilde{c}*(\tau^2)^(1/2)/\theta-log(exp((\tau^2)^(1/2)*\tilde{c}/\theta)+exp(\tilde{c})))}
#' 			\eqn{\tilde{c}} is chosen such that \eqn{P(\tau^2<=\theta)>=0.95}.
#'
#' @references Nadja Klein and Thomas Kneib (2015). Scale-Dependent Priors for Variance Parameters in Structured Additive Distributional Regression. 
#' \emph{Working Paper}.
#'
#' @seealso \code{\link{rapprox_unif}},\code{\link{dapprox_unif}}
#' @export
 
 
papprox_unif <- function(x,scale,tildec=13.86294) 
  {
  arg <- x^(0.5)*tildec/scale 
  #ret <- (sqrt(x)*tildec - log(1+exp(arg-tildec)) + log(1+exp(-tildec))) / (scale*(1+log(exp(-tildec))+tildec))
  ret <- ( tildec*x^(0.5)/scale - log(exp(arg)+exp(tildec))) / ((tildec + log(1+exp(-tildec)) ))+1
  return(ret)
  }
  

NULL

dapprox_unif_exp <- function(x,scale,tildec=13.86294,log=FALSE) 
  {
  ex <- exp(x)
  arg <- ex^(0.5)*tildec/scale - tildec
  if(log)
    ret <- 0.5*x-log(2)-log(scale)-log(1+log(1+exp(-tildec))/tildec)-log(1+exp(arg))
  else
    ret <- 0.5*sqrt(ex)*(1/(1+exp(arg)))/(scale*(1+log(1+exp(-tildec))/tildec))
  return(ret)
  }  
  
  
NULL

#' Computing Designmatrix for Splines
#' 
#' This function computes the design matrix for Bayesian P-splines as it would be done in 
#' BayesX. The implementation currently on works properly for default values (knots=20, degree=3).

#' 
#' @param x the covariate vector.  
#' @param degree of the B-splines, default is 3.
#' @param knots number of knots, default is 20. 
#' @param min_x the left interval boundary, default is min(x). 
#' @param max_x the right interval boundary, defalut is max(x). 
#' @return a list with design matrix at distinct covariates, design matrix at all observations, 
#'  index of sorted observations, the difference matrix, precision matrix and the knots used.
#' @author Nadja Klein
#' @references Stefan Lang and Andy Brezger (2004). Bayesian P-Splines.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{13}, 183--212.
#' 
#' Belitz, C., Brezger, A., Klein, N., Kneib, T., Lang, S., Umlauf, N. (2015): BayesX - Software for Bayesian inference in structured additive regression models.
#' Version 3.0.1. Available from http://www.bayesx.org. 
#' @import splines
#' @export

DesignM <- function(x,degree=3,knots=20,min_x=min(x),max_x=max(x)) {

  ##############################################################################################
  #reduce establishment of design matrix to unique and sorted observations
  x_work <- unique(sort.int(x,index.return=TRUE)$x)

  min_x_work <- min_x
  max_x_work <- max_x
  #find indices to get back to whole covariate vector of length n(number of observations.)
  ix <- c()
  for(i in 1:length(x)) {
    ix <- c(ix,which(x_work==x[i]))
  }

  ##############################################################################################  
  #set degree of B-splines and number of inner knots
  degree <- degree
  knots <- knots
  l <- degree
  m <- knots
  hilfsgr <- (max_x_work-min_x_work)/100
  step <- (max_x_work+hilfsgr-min_x_work+hilfsgr)/(knots-1)
  knoten_B=seq(min_x_work-hilfsgr-degree*step,max_x_work+hilfsgr+degree*step,by=step)
  if(length(knoten_B)>1) {
    Z_BayesX<- splineDesign(knoten_B,x_work,ord=degree+1,outer.ok=T) #designmatrix of x_work as in BayesX
  } else {
	knoten=seq(min_x_work-0.001-(degree)*((max_x_work+0.001-min_x_work+0.001)/(knots-1)),max_x_work+0.001+(degree)*((max_x_work+0.002-min_x_work)/(knots-1)),length=26)
	Z_BayesX<- splineDesign(knoten,x_work,ord=degree+1,outer.ok=T) #designmatrix of x_work as in BayesX
  }
 
  Zcomp_BayesX <- Z_BayesX[ix,]
  #Construction of matrix of second differences (beta_j-2*beta_{j-1}+beta_{j-2})
  Dmatrix <- matrix(0,nrow=m+l-3,ncol=m+l-1) 
  for(i in 1:(m+l-3)){
    Dmatrix[i,i] <- 1
    Dmatrix[i,i+1] <- -2
    Dmatrix[i,i+2] <- 1
  }
  
  #Kmatrix, product of the transposed Dmatrix and the Dmatrix (needed for penalizing term)
  Kmatrix <- t(Dmatrix)%*%Dmatrix 
  
  return(list(Z_B=Z_BayesX,Zcomp_B=Zcomp_BayesX,ix=ix,Dmatrix=Dmatrix,Kmatrix=Kmatrix,knoten_B=knoten_B))
  
  
}