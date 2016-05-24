#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       give.init.y.R                       ####
####                                                 ####
#### FUNCTIONS:  give.init.y                         ####
#########################################################

### ======================================
### give.init.y
### ======================================
###
### Compute initial augmented observations
### * used by functions that exploit G-splines
###
### 15/01/2005: 'give.init.y'
### 19/01/2005: 'give.init.y2'
### ===========================================
###
### \item{init.y}{initial (augmented) observations possibly given by the user.                        
###       They are partially checked for consistency and these supplied by the user                   
###       used in the resulting object. This should be either vector of length \eqn{n}{n}             
###       where \eqn{n}{n} is a~sample size if the dimension is one or a~matrix with                  
###       2 columns and \eqn{n}{n} rows if the dimension is two.}                                     
### \item{dim}{dimension of the G-spline, 1 or 2.}                                                    
### \item{y.left}{observed, left or right censored log(event time) or the lower limit                 
###       of the interval censored observation. Sorted in a~transposed order compared                 
###       to \code{init.y}.}                                                                          
### \item{y.right}{upper limit of the interval censored observation, whatever if the observation      
###       is not interval-censored sorted in a~transposed order compared to \code{init.y}.}           
### \item{status}{status indicator vector/matrix. 1 for observed times, 0 for right censored times,   
###       2 for left censored times, 3 for interval censored times.}                                  
###
### \value{
###   a~vector or matrix with the same structure as \code{init.y}, i.e. with 2~columns and
###   \eqn{n}{n} rows in the case of the bivariate data.
### }
give.init.y <- function(init.y, dim, y.left, y.right, status)
{
  n <- length(status)/dim
  
  if (missing(init.y)) init.y <- NULL
  if (!length(init.y)) init.y <- rep(NA, n*dim)
  else{
    if (dim == 1){
      if (length(init.y) != n)  stop("init$y has a different length than the data")
    }
    else{
      if (is.null(dim(init.y))) stop("init$y has dimensions inconsistent with the data")
      else                      if (dim(init.y)[1] != n | dim(init.y)[2] != 2) stop("init$y has dimensions inconsistent with the data")
    }                          
  }    
  ysm <- as.vector(t(init.y))     ## transposition to be consistent with y.left and y.right    

  not.na <- !is.na(ysm)
  if (sum(ysm[not.na & status == 0] < y.left[not.na & status == 0])) stop("Inconsistent init$y found (lower than right censored obs.)")
  if (sum(ysm[not.na & status == 2] > y.left[not.na & status == 2])) stop("Inconsistent init$y found (higher than left censored obs.)")
  if (sum(ysm[not.na & status == 3] < y.left[not.na & status == 3])) stop("Inconsistent init$y found (outside the interval)")
  if (sum(ysm[not.na & status == 3] > y.right[not.na & status == 3])) stop("Inconsistent init$y found (outside the interval)")      
  ysm[status == 1] <- y.left[status == 1]
  ysm[!not.na & status == 0] <- y.left[!not.na & status == 0]
  ysm[!not.na & status == 2] <- y.left[!not.na & status == 2]
  ysm[!not.na & status == 3] <- 0.5*(y.left[!not.na & status == 3] + y.right[!not.na & status == 3])
  init.y <- matrix(ysm, ncol=dim, nrow=n, byrow=TRUE)
  
  return(init.y)
}

###############################################################################################################

### Version for bayesBisurvreg
### * mid-point imputation performed using C++ functions
### =====================================================
give.init.y2 <- function(init.y, init2.y, dim, design, design2, doubly)
{
  thispackage <- "bayesSurv"
  #thispackage <- NULL
  
  t1.left <- matrix(design$Y[, 1], nrow=dim)
  if (design$nY > 2) t1.right <- matrix(design$Y[, 2], nrow=dim)
  else               t1.right <- matrix(rep(1, length(t1.left)), nrow=dim)
  status1 <- matrix(design$Y[, design$nY], nrow=dim)
  nP <- length(t1.left)
  n <- nP/dim
  
  ### DOUBLY censored observations
  ### ----------------------------
  if (doubly){
    t2.left <- matrix(design2$Y[, 1], nrow=dim)
    if (design2$nY > 2) t2.right <- matrix(design2$Y[, 2], nrow=dim)
    else                t2.right <- matrix(rep(1, length(t2.left)), nrow=dim)
    status2 <- matrix(design2$Y[, design2$nY], nrow=dim)

    midImp <- .C("midimputeDataDoubly", err=integer(1), t1=double(nP), t2=double(nP), as.integer(nP),
                                        as.double(t1.left), as.double(t1.right), as.integer(status1),
                                        as.double(t2.left), as.double(t2.right), as.integer(status2), PACKAGE = thispackage)
    if (midImp$err) stop("Inconsistent data supplied")
    t1 <- midImp$t1
    t2 <- midImp$t2

      ## Manipulation with init2.y
    if (missing(init2.y)) init2.y <- rep(NA, nP)
    if (!length(init2.y)) init2.y <- rep(NA, nP)
    else{
      if (dim == 1){
        if (length(init2.y) != n)  stop("init2$y has a different length than the data")
      }
      else{
        if (is.null(dim(init2.y))) stop("init2$y has dimensions inconsistent with the data")
        else                       if (dim(init2.y)[1] != n | dim(init2.y)[2] != 2) stop("init2$y has dimensions inconsistent with the data")
      }                          
    }    
    ysm2 <- as.vector(t(init2.y))                           ## transposition to be consistent with t2.left, t2.right and t2
    ysm2[is.na(ysm2)] <- log(t2[is.na(ysm2)])               ## logarithmic transformation for non-user supplied values
    init2.y <- matrix(ysm2, ncol=dim, nrow=n, byrow=TRUE)
  }

  ### SIMPLE censoring
  ### ----------------
  else{
    midImp <- .C("midimputeData", err=integer(1), t1=double(nP), as.integer(nP),
                                  as.double(t1.left), as.double(t1.right), as.integer(status1), PACKAGE = thispackage)
    if (midImp$err) stop("Inconsistent data supplied")
    t1 <- midImp$t1

    t2.left <- 0
    t2.right <- 0
    status2 <- 0
    init2.y <- 0
  }    

    ## Manipulation with init.y
  if (missing(init.y)) init.y <- rep(NA, nP)
  if (!length(init.y)) init.y <- rep(NA, nP)
  else{
    if (dim == 1){
      if (length(init.y) != n)  stop("init$y has a different length than the data")
    }
    else{
      if (is.null(dim(init.y))) stop("init$y has dimensions inconsistent with the data")
      else                      if (dim(init.y)[1] != n | dim(init.y)[2] != 2) stop("init$y has dimensions inconsistent with the data")
    }                          
  }    
  ysm <- as.vector(t(init.y))                         ## transposition to be consistent with t1.left, t1.right and t1
  ysm[is.na(ysm)] <- log(t1[is.na(ysm)])              ## logarithmic transformation for non-user supplied values
  init.y <- matrix(ysm, ncol=dim, nrow=n, byrow=TRUE)

  y1.left <- log(t1.left)
  y1.right <- log(t1.right)

  toreturn <- list(init.y=init.y, init2.y=init2.y, y1.left=y1.left, y1.right=y1.right, status1=status1,
                                                   t2.left=t2.left, t2.right=t2.right, status2=status2)

  return(toreturn)  
}  
