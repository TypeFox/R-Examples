############################################################################################################################################
##################################               S U M M A R Y    S T A T I S T I C S                     ##################################
############################################################################################################################################
############ All individual and grouped summary statistics are stored here


##########################
######## INDIVIDUAL SUMMARY STATISTICS
##########################
#' Estimate auto-covariances for multiple datasets.
#' 
#' @description Function that, give time series data, transforms them into auto-covariances with different lags.
#'              
#' @param x a matrix. Each column contains a replicate series.
#' @param max.lag How many lags to use.    
#' 
#' @return a matrix where each column contains the coefficients for a different replicate. The first coefficient
#'         corresponds to lag == 0, hence it is the variance, the second is the covariance one step ahead and so 
#'         on.
#' @author Simon N. Wood, maintainer Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @rdname slAcf
#' @export
#' @examples
#' library(synlik)
#' set.seed(10)
#' x <- matrix(runif(1000),100,10)
#' acf <- slAcf(x)

slAcf <- function(x, max.lag=10) {
  ## `x' is a matrix containing replicate simulations in its columns.
  ## slAcf turns these into acf's
  
  NAcode <- -1e70
  x[is.na(x)] <- NAcode
  
  acf <- matrix(0,max.lag+1,ncol(x))
  oo <- .C("slacf",acf=as.double(acf),x=as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
           as.integer(max.lag),as.double(NAcode),correlation=as.integer(0),PACKAGE="synlik")
  
  acf <- matrix(oo$acf,max.lag+1,ncol(x))
  acf[acf == NAcode] <- NA
  acf
  
} 

#' Estimate non-linear autoregressive coefficients
#' 
#' @description Function that, give time series data, transforms them into summary statistics
#'              using polynomial autoregression.
#'              
#' @param x a matrix. Each column contains a replicate series.
#' @param lag vector of lags, for rhs terms.
#' @param power vector of powers, for rhs terms.       
#' 
#' @return a matrix where each column contains the coefficients for a different replicate.
#' @author Simon N. Wood, maintainer Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @rdname nlar
#' @export
#' @examples
#'   library(synlik)
#'   set.seed(10)
#'   x <- matrix(runif(200),100,2)
#'   beta <- nlar(x,lag=c(1,1),power=c(1,2))
#'   y <- x[,1]
#'   y <- y - mean(y)
#'   z <- y[1:99];y <- y[2:100]
#'   lm(y~z+I(z^2)-1)
#'   beta
#'   
#'   ## NA testing
#'   x[5,1] <- x[45,2] <- NA
#'   beta <- nlar(x,lag=c(1,1),power=c(1,2))
#'   y <- x[,1]
#'   y <- y - mean(y,na.rm=TRUE)
#'   z <- y[1:99];y <- y[2:100]
#'   lm(y~z+I(z^2)-1)
#'   beta
#'   
#'   ## higher order...
#'   set.seed(10)
#'   x <- matrix(runif(100),100,2)
#'   beta <- nlar(x,lag=c(6,6,6,1,1),power=c(1,2,3,1,2))
#'   k <- 2
#'   y <- x[,k]
#'   y <- y - mean(y)
#'   ind <- (1+6):100
#'   y6 <- y[ind-6];y1 <- y[ind-1];y <- y[ind]
#'   beta0 <- coef(lm(y~y6+I(y6^2)+I(y6^3)+y1+I(y1^2)-1))
#'   as.numeric(beta[,k]);beta0;beta0-as.numeric(beta[,k])
  
nlar <- function(x,lag,power) {
  ## relatively efficient polynomial autoregression for multiple reps.
  ## each column of `x' is a replicate. 
  ## `lag[i]' is the lag for term i on rhs of autoregression
  ## `power[i]' is the power for term i on rhs of autoregression 
  beta <- matrix(0,length(lag),ncol(x))
  
  NAcode <- -1e70
  x[is.na(x)] <- NAcode  
  
  oo <- .C("slnlar",beta = as.double(beta), x = as.double(x),
           n=as.integer(nrow(x)),n.reps=as.integer(ncol(x)),n.terms=as.integer(length(lag)),
           as.integer(lag),as.integer(power),as.double(NAcode),PACKAGE="synlik")
  
  beta <- matrix(oo$beta,length(lag),ncol(x))
  
  beta
} 

#' Summarize marginal distribution of (differenced) series.
#' 
#' @description Summarizes (difference) distribution of replicate series, by regressing ordered 
#'              differenced series on a reference series (which might correspond to observed data).
#'              
#' @param x a matrix. Each column contains a replicate series.
#' @param z vector of lags, for rhs terms.    
#' @param np maximum power on rhs of regression.
#' @param diff order of differencing (zero for none).
#' 
#' @return a matrix where each column contains the coefficients for a different replicate.
#' @author Simon N. Wood, maintainer Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @rdname orderDist
#' @export
#' @examples
#' library(synlik)
#' set.seed(10)
#' n <- 100;nr <- 3
#' x <- matrix(runif(n*nr),n,nr)
#' z <- runif(n)
#' beta <- orderDist(x,z,np=3,diff=1)
#' 
#' zd <- z;xd <- x[,3]
#' zd <- diff(zd,1);xd <- diff(xd,1)
#' zd <- sort(zd);zd <- zd - mean(zd)
#' xd <- sort(xd);xd <- xd - mean(xd)
#' lm(xd~zd+I(zd^2)+I(zd^3)-1)

orderDist <- function(x,z,np=3,diff=1) {
  ## Routine to obtain coefficients summarizing distribution of (differenced) columns
  ## of x, by regression of sorted differenced columns of x on sorted differenced z's. 
  ## regression is with order `np' polynomial (no intercept as all centred). `diff'
  ## is order of differencing to apply.
  
  beta <- matrix(0,np,ncol(x))
  oo <- .C("order_reg",beta=as.double(beta), as.double(x),as.double(z),n=as.integer(nrow(x)),
           as.integer(ncol(x)),as.integer(np),as.integer(diff),PACKAGE="synlik")
  
  beta <- matrix(oo$beta,np,ncol(x))
  
} 




##########################
######## Grouped SUMMARY STATISTICS
##########################

########
### WOOD2010 original statistics
########

rickerStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(orderDist(tx, obsData, np=3, diff=1))        ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(tx^.3, lag=c(1,1), power=c(1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowSums(x==0))         ## mean values of Y, # of 0's
  X0 <- cbind(X0, t(slAcf(tx, max.lag=5)))                 ## autocovariances up to lag 5 (the first element is the variance)
  
  X0
}


#########################
####### Blowfly Statistics
#########################

blowStats <- function(x, extraArgs, ...){
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  tx <- t(x)
  
  X0 <- t(orderDist(t(x), obsData,np=3,diff=1)) ## cubic regs coeff of ordered diffs
  X0 <- cbind(X0, t(nlar(t(x),lag=c(6, 6, 6, 1, 1),power=c(1,2,3,1,2)))) ## least square coeff
  X0 <- cbind(X0, rowMeans(x), rowMeans(x) - apply(x,1,median)) 
  X0 <- cbind(X0, t(slAcf(t(x), max.lag=11)))  ## autocovariances up to lag 5 (the first element is the variance)
  
  X0 <- cbind(X0, apply( t( abs( apply( sign( apply(x,1,diff) ), 2, diff ) )), 1, sum)/2) #Number of turning points 
  
  X0
}