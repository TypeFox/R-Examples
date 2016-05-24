#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       give.init.r.R                       ####
####                                                 ####
#### FUNCTIONS:  give.init.r                         ####
#########################################################

### ======================================
### give.init.r
### ======================================
###
### Compute initial allocations
### * used by functions that exploit G-splines
###
### 15/01/2005
### ===========================================
###
### \item{init.r}{initial allocations possibly given by the user. This should be a~vector of length                       
###    \eqn{n}{n} where n is a~sample size or a matrix with 2 columns if dim=2.}
### \item{init.y}{correctly computed initial observations the G-spline is estimated from. In the case of regression       
###    these should be replaced by residuals. This must be either a~vector or matrix (in the format as returned           
###    by the function \code{\link{give.init.y}}).}                                                                       
### \item{dim}{dimension of the G-spline, 1 or 2.}                                                                        
### \item{KK}{vector of length \code{dim} with \eqn{K}{K} coefficients defining the G-spline.}                            
### \item{gamma}{vector of length \code{dim} with initial \eqn{\gamma}{gamma} parameters of the G-spline.}                
### \item{sigma}{vector of length \code{dim} with initial \eqn{\sigma}{sigma} parameters of the G-spline.}                
### \item{c4delta}{vector of length \code{dim} with constants to compute the distance between two knots                   
###    defining the G-spline.}                                                                                            
### \item{intcpt}{vector of length \code{dim} with initial values of the intercept term of the G-spline.}                 
### \item{scale}{vector of length \code{dim} with initial values of the scale parameters of the G-spline.}                
###
give.init.r <- function(init.r, init.y, dim, KK, gamma, sigma, c4delta, intcpt, scale)
{
  thispackage <- "bayesSurv"
  #thispackage <- NULL

  n <- length(init.y)/dim
  if (dim == 1) init.y <- matrix(init.y, ncol=1)

  if (missing(init.r)) init.r <- NULL
  if(!length(init.r)){  
    for (j in 1:dim){
      knots <- gamma[j] + seq(-KK[j], KK[j], by=1)*c4delta[j]*sigma[j]
      ind <- .C("findClosestKnot", ind=integer(n), as.double(knots), as.double((init.y[,j] - intcpt[j])/scale[j]),
                                   as.integer(length(knots)), as.integer(n),
                                   PACKAGE = thispackage)
      ind <- ind$ind
      if (j == 1) init.r <- ind - KK[j]
      else        init.r <- cbind(init.r, ind - KK[j])
    }    
  }
  if (dim == 1){
    if (length(init.r) != n)  stop("init$r has a different length than the data")
  }
  else{
    if (is.null(dim(init.r))){ stop("init$r has dimensions inconsistent with the data") }
    else                     { if (dim(init.r)[1] != n | dim(init.r)[2] != 2) stop("init$r has dimensions inconsistent with the data") }
  }                        
  if (sum(is.na(init.r))) stop("Missing values found in init$r")    
  if (dim == 1){ if (sum(init.r < -KK[1]) | sum(init.r > KK[1])) stop("init$r inconsistent with prior$K") }
  else{
    if (sum(init.r[,1] < -KK[1]) | sum(init.r[,1] > KK[1])) stop("init$r inconsistent with prior$K[1]")
    if (sum(init.r[,2] < -KK[2]) | sum(init.r[,2] > KK[2])) stop("init$r inconsistent with prior$K[2]")
  }

  return(init.r)  
}
