##########################################################################
##             Archimedean Liouville generation and estimation          ##
##                        Version of August 16th, 2015                  ##
##                     Last modified - August 16th, 2015                ##
##                         Built - February 18th, 2015                  ##
##########################################################################

# Code adapted from GumbelLiouville.R (Neslehova and Genest - 2013)
# and from LiouvilleFunction.R (McNeil and Neslehova - 2010)
# using functions implemented by Hofert, Machler and McNeil (2012)
# in the package copula

# Warning: different definition used by latter for the Clayton
# The implementation below uses the generator definition found in QRM

# The Archimedean families implemented are Clayton, Gumbel, Frank, 
# Ali-Mikhail-Haq (abbreviated AMH) and Joe
# require(copula)

##########################################################################


## Necessary global imports for the methods for NAMESPACE
## Note that stats cannot be imported because it is already loaded by copula
#' @useDynLib lcopula
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom graphics "axis" "lines" "mtext" "plot.window" "points"
#' @importFrom stats "integrate" "optimise" "pbeta" "quantile" "rgamma" "runif" "uniroot"
#' @import copula
#' @importFrom utils combn
#' @importFrom pcaPP cor.fk
#' @importFrom Rcpp evalCpp
NULL


#' Archimedean copula sampler
#'
#' Sampler based on the Marshall-Olkin algorithm
#'
#' @param n sample size
#' @param family family of the Archimedean copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param d dimension of sample
#' @param theta parameter of the Archimedean copula
#' @return a sample of dimension \code{n} by \code{d} from the Archimedean copula
#' @export
#' @examples
#' #Sample from a Gumbel Archimedean copula
#' rarchi(n = 100, "gumbel", d = 4, theta = 2)
#' #Sample from the independence copula
#' rarchi(n = 100, "gumbel", d = 4, theta = 1)
rarchi <-  function(n, family, d, theta){
  family <-  match.arg(family, c("clayton","gumbel","frank","AMH","joe"))
  #Marshall-Olkin algorithm
  illegalpar <-  switch(family, 
                       clayton = copClayton@paraConstr(theta), 
                       gumbel = copGumbel@paraConstr(theta), 
                       frank = copFrank@paraConstr(theta), 
                       AMH = copAMH@paraConstr(theta), 
                       joe = copJoe@paraConstr(theta))
  if(!illegalpar){
    stop("Illegal parameter value")
  }
  independence <-  switch(family, 
                         clayton = (theta == 0), 
                         gumbel = (theta == 1), 
                         frank = (theta == 0), 
                         AMH = (theta == 0), 
                         joe = (theta == 1))
  U <-  runif(n * d)
  U <-  matrix(U, nrow = n, ncol = d)
  if(independence)
    return(U)
  Y <-  switch(family, 
              clayton = rgamma(n, shape = 1/theta, scale = theta), 
              gumbel = copGumbel@V0(n, theta), 
              frank = copFrank@V0(n, theta), 
              AMH = copAMH@V0(n, theta), 
              joe = copJoe@V0(n, theta))
  Y <-  matrix(Y, nrow = n, ncol = d)
  phi.inverse <-  switch(family, 
                        clayton = function(t, theta){(1 + theta*t)^(-1/theta)}, 
                        gumbel = copGumbel@psi, 
                        frank = copFrank@psi, 
                        AMH = copAMH@psi, 
                        joe = copJoe@psi)
  return(phi.inverse( - log(U)/Y, theta))
}


##########################################################################
#' Liouville copula sampler
#'
#' @param n sample size
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers) Specifies (implictly) the dimension of sample
#' @param theta parameter of the corresponding Archimedean copula
#' @param reverse if \code{TRUE}, return sample from the corresponding survival copula
#' @return a sample of dimension \code{n} by \code{length(alphavec)} from the Liouville copula
#' @export
#' @examples
#' # Sample from a Gumbel Liouville copula
#' rliouv(n = 100, family = "gumbel", alphavec = c(2, 3), theta = 2)
rliouv <-  function(n = 100, family, alphavec, theta, reverse = FALSE){
  family <-  match.arg(family, c("clayton","gumbel","frank","AMH","joe"))
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  alpha <-  sum(alphavec)
  d <-  length(alphavec)
  illegalpar <-  switch(family, 
                       clayton = copClayton@paraConstr(theta), 
                       gumbel = copGumbel@paraConstr(theta), 
                       frank = copFrank@paraConstr(theta), 
                       AMH = copAMH@paraConstr(theta), 
                       joe = copJoe@paraConstr(theta))
  if(!illegalpar)
    stop("Illegal parameter value")
  t <-  rarchi(n, family, alpha, theta)
  Xtildedata <-  switch(family, 
                       clayton = (t^(-theta)-1)/theta, 
                       gumbel = copGumbel@iPsi(t, theta), 
                       frank = copFrank@iPsi(t, theta), 
                       AMH = copAMH@iPsi(t, theta), 
                       joe = copJoe@iPsi(t, theta))
  groups <-  c(0, cumsum(alphavec))
  Xdata <-  matrix(0, nrow = n, ncol = d)
  Udata <-  Xdata
  for (j in 1:d){
    class <-  ((groups[j]+1):groups[j+1])
    for (i in class)
      Xdata[, j] <-  Xdata[, j]+Xtildedata[, i]
  }
  for (j in 1:d){
    Udata[, j] <-  sliouvm(Xdata[, j], family, alphavec[j], theta)
  }
  if(reverse == F){
    return(Udata)
  } else{
    return(1-Udata)
  }
}

##########################################################################

#' Marginal survival function for Liouville copulas
#'
#'Default behavior as in original GumbelLiouville.R
#'
#' @param x sample from copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alpha marginal allocation parameter (must be integer)
#' @param theta parameter of the corresponding Archimedean copula
#' @return a vector of same length as \code{x} with the survival probabilities
#' @export
#' @examples 
#' x <- rliouv(n = 100, family = "gumbel", alphavec <- c(2,3), theta = 2)
#' sliouvm(x[,1], family="gumbel", alpha=alphavec[1], theta=2)
sliouvm <- function(x, family, alpha, theta){
  family <-  match.arg(family, c("clayton","gumbel","frank","AMH","joe"))
  alpha <- as.integer(alpha)
  if(alpha==0){stop("Invalid parameter alpha")}
  out <- switch(family, 
               gumbel = copGumbel@psi(t = x, theta), 
               clayton = copClayton@psi(t = x*theta, theta), 
               frank = copFrank@psi(t = x, theta), 
               AMH = copAMH@psi(t = x, theta), 
               joe = copJoe@psi(t = x, theta))
  if(alpha > 1){
    for(i in 1:(alpha-1)){
      out <-  out+ exp(i*log(x)-lgamma(i+1)+
                        switch(family, 
                               gumbel = copGumbel@absdPsi(x, theta = theta, degree = i, log = T), 
                               clayton = copClayton@absdPsi(x*theta, theta = theta, degree = i, log = T)+i*log(theta), 
                               frank = copFrank@absdPsi(x, theta = theta, degree = i, log = T), 
                               AMH = copAMH@absdPsi(x, theta = theta, degree = i, log = T), 
                               joe = copJoe@absdPsi(x, theta = theta, degree = i, log = T))
      )
    }
  }
  return(out)
}

##########################################################################
#' Marginal distribution function for Liouville copula
#'
#' Calls \code{sliouvm} to return cdf
#'
#' @param x sample from copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alpha marginal allocation parameter (must be integer)
#' @param theta parameter of the corresponding Archimedean copula
#' @export
#' @examples 
#' x <- rliouv(n = 100, family = "gumbel", alphavec <- c(2,3), theta = 2)
#' pliouvm(x[,1], family="gumbel", alpha=alphavec[1], theta=2)
#' 
#' @return a vector of same length as \code{x} with the quantiles
pliouvm <- function(x, family, alpha, theta){
  family <-  match.arg(family, c("clayton","gumbel","frank","AMH","joe"))
  return(1-sliouvm(x, family, alpha, theta))
}

##########################################################################
#' Multiple marginal survival function for Liouville copulas
#'
#' To apply different parameters (i.e. treat margins simultaneously)
#'
#' @param x sample from copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param theta parameter of the corresponding Archimedean copula
#' @export
#' @examples
#' x <- rliouv(n = 100, family = "gumbel", alphavec <- c(2,3), theta = 2)
#' sliouv_m(x, family="gumbel", alphavec=c(2,3), theta=2)
#' all(sliouv_m(x, family="gumbel", alphavec=c(2,3), theta=2)[,1]-
#'   sliouvm(x[,1], family="gumbel", alpha=2, theta=2)==0)
#' @return a matrix of same length as \code{x} with the survival probabilities
sliouv_m <- function(x, family, alphavec, theta){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton","gumbel","frank","AMH","joe"))
  eval.Surv <-  function(x, alpha){
    der <-  matrix(nrow = alpha, ncol = length(x))
    der[1, ] <-  switch(family, 
                      gumbel = copGumbel@psi(t = x, theta = theta), 
                      clayton = copClayton@psi(t = x*theta, theta = theta), 
                      frank = copFrank@psi(t = x, theta = theta), 
                      AMH = copAMH@psi(t = x, theta = theta), 
                      joe = copJoe@psi(t = x, theta = theta))
    if(alpha > 1){
      for(i in 1:(alpha-1)){
        der[i+1, ] <-  exp(i*log(x)-lgamma(i+1)+switch(family, 
                                                     gumbel = copGumbel@absdPsi(x, theta = theta, degree = i, log = T), 
                                                     clayton = copClayton@absdPsi(x*theta, theta = theta, degree = i, log = T)+i*log(theta), 
                                                     frank = copFrank@absdPsi(x, theta = theta, degree = i, log = T), 
                                                     AMH = copAMH@absdPsi(x, theta = theta, degree = i, log = T), 
                                                     joe = copJoe@absdPsi(x, theta = theta, degree = i, log = T)))
      }
    }
    return(colSums(der))
  }
  if(is.matrix(x)){
    if(dim(x)[2] != length(alphavec))
      stop("Length of the parameter vector does not match input")
    result <- list()
    for(j in 1:(dim(x)[2])){
      result[[j]] <- eval.Surv(x[, j], alphavec[j])
    }
    result <- unlist(result)
  }  else{
    result <-  eval.Surv(x, alphavec[1])
  }
  return(matrix(result, nrow = nrow(x), ncol = ncol(x)))
}

##########################################################################
#' Marginal inverse survival function for Liouville copulas
#'
#' Survival analog to \code{'q'} in R for functions
#'
#' @param u vector of survival probabilities
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alpha marginal allocation parameter (must be integer)
#' @param theta parameter of the corresponding Archimedean copula
#' @export
#' @return a vector of same length as \code{u} with the quantile at \code{1-u}
#' @examples
#' u <- rliouv(n = 100, family = "clayton", alphavec <- c(2,3), theta = 2)
#' isliouvm(u=u[,1], family="clayton", alpha=2, theta=2)
isliouvm <-  function(u, family, alpha, theta){
  alpha <- as.integer(alpha)
  if(alpha==0){stop("Invalid parameter alpha")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  out <-  rep(0, length(u))
  rootfunc <-  function(z, u, family, theta, alpha)
  {
    sliouvm(z, family, alpha, theta)-u
  }
  for (i in 1:length(u))
  {
    #Handling of 0 for the Gumbel
    tmp <-  uniroot(rootfunc, lower = 0, upper = 1e+60, u = u[i], f.lower = 1, f.upper = -1e+60, 
                   family = family, theta = theta, alpha = alpha)
    out[i] <- tmp$root
  }
  out
}

##########################################################################
#' Multiple marginal inverse survival function for Liouville copulas
#'
#' @param u vector of survival probabilities
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param theta parameter of the corresponding Archimedean copula
#' @return a vector of same length as \code{u} with the quantile at 1-u
#' @export
#' @examples
#' u <- rliouv(n = 100, family = "clayton", alphavec <- c(2,3), theta = 2)
#' isliouvm_m(u=u, family="clayton", alphavec=c(2,3), theta=2)
isliouvm_m <-  function(u, family, alphavec, theta){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  out <-  matrix(nrow = nrow(u), ncol = ncol(u))
  rootfunc <-  function(z, u, family, theta, alpha)
  {
    sliouvm(z, family, alpha, theta)-u
  }
  for (i in 1:nrow(u)){
    for(j in 1:ncol(u)){
      tmp <-  uniroot(rootfunc, lower = 0, upper = 1e+60, u = u[i, j], f.lower = 1, f.upper = -1e+60, 
                     family = family, theta = theta, alpha = alphavec[j])
      out[i, j] <- tmp$root
    }
  }
  out
}

##########################################################################
#' Joint survival function for Liouville copulas
#'
#' @param theta parameter of the corresponding Archimedean copula
#' @param data sample matrix from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#'
#' @return a matrix of survival probabilities
#' @export
#' @examples
#' x <- rliouv(n = 100, family = "clayton", alphavec <- c(2,3), theta = 2)
#' sliouv(data=x, family="clayton", alphavec=c(2,3), theta=2)
sliouv <- function(theta, data, family, alphavec){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  combo = unique(combn(unlist(sapply(alphavec, function(y){seq(1:y)})), length(alphavec)), MARGIN = 2)-1
  result = apply(combo, 2, function(alpha_combo){
    apply(data, 1, function(x){
      ifelse(sum(alpha_combo) != 0, exp(sum(alpha_combo*log(x))-sum(lgamma(alpha_combo+1))+
                                       switch(family, 
                                              gumbel = copGumbel@absdPsi(sum(x), theta = theta, degree = sum(alpha_combo), log = T), 
                                              clayton = copClayton@absdPsi(t = theta*sum(x), theta = theta, degree = sum(alpha_combo), log = T)+
                                                sum(alpha_combo)*log(theta), 
                                              frank = copFrank@absdPsi(sum(x), theta = theta, degree = sum(alpha_combo), log = T), 
                                              AMH = copAMH@absdPsi(sum(x), theta = theta, degree = sum(alpha_combo), log = T), 
                                              joe = copJoe@absdPsi(sum(x), theta = theta, degree = sum(alpha_combo), log = T))), 
             switch(family, 
                    gumbel = copGumbel@psi(sum(x), theta = theta), 
                    clayton = copClayton@psi(t = theta*sum(x), theta = theta), 
                    frank = copFrank@psi(sum(x), theta = theta), 
                    AMH = copAMH@psi(sum(x), theta = theta), 
                    joe = copJoe@psi(sum(x), theta = theta)))
    }
    )
  }
  )
  return(rowSums(result))
}

##########################################################################
#' Copula function for Liouville copulas
#'
#' @param pseudo sample matrix from a Liouville copula
#' @param theta parameter of the corresponding Archimedean copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#'
#' @return a matrix of survival probabilities
#' @export
#' @examples
#' x <- rliouv(n=100, family="frank", theta=1.5, alphavec=c(2,3))
#' pliouv(theta=1.5, pseudo=x,family="frank", alphavec=c(2,3))
pliouv <- function(pseudo, theta, family, alphavec){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  sliouv(theta, data = isliouvm_m(pseudo, family, alphavec, theta), family, alphavec)
}

##########################################################################
#' Copula function for Liouville copulas
#'
#' @param x sample matrix from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param theta parameter of the corresponding Archimedean copula
#' @param is.log if \code{TRUE}, will return the log-likelihood value
#'
#' @return value of multivariate density
#' @export
#' @examples
#' x <- rliouv(n = 100, family = "clayton", alphavec <- c(2,3), theta = 2)
#' dliouv(x=x, family="clayton", alphavec=c(2,3), theta=2, TRUE)
dliouv <- function(x, family, alphavec, theta, is.log = F){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  if(!dim(x)[2] == length(alphavec))
  {
    stop("Invalid argument: alpha and x are not of same length")
  }
  #Check for validity of the returned parameter in light of family
  alpha <-  sum(alphavec)
  x.norm <-  rowSums(x)
  psi.alpha <-  switch(family, 
                      gumbel = copGumbel@absdPsi(x.norm, theta = theta, degree = alpha, log = T), 
                      clayton = copClayton@absdPsi(x.norm*theta, theta = theta, degree = alpha, log = T)+alpha*log(theta), 
                      frank = copFrank@absdPsi(x.norm, theta = theta, degree = alpha, log = T), 
                      AMH = copAMH@absdPsi(x.norm, theta = theta, degree = alpha, log = T), 
                      joe = copJoe@absdPsi(x.norm, theta = theta, degree = alpha, log = T))
  arg1 <-  t(apply(x, 1, function(x){(alphavec-1)*log(x)}))
  return_val <-  sum(arg1)+sum(psi.alpha)-dim(x)[1]*sum(lgamma(alphavec))
  if(is.log == T){return(return_val)
  } else{
    return(exp(return_val))
  }
}


##########################################################################
#' Marginal density function for Liouville copulas
#'
#'The function is using the telescoping series property - that is, it computes only the last derivative.
#'
#' @param x sample vector from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alpha allocation parameter (must be an integer)
#' @param theta parameter of the corresponding Archimedean copula
#'
#' @return value of marginal density
#' @export
#' @examples
#' samp <- rliouv(n = 100, family = "clayton", alphavec <- c(2,3), theta = 2)
#' dliouvm(x=samp[,1], family="clayton", alpha=2, theta=2)
#' sum(log(dliouvm(x=samp[,1], family="clayton", alpha=2, theta=2)))
dliouvm <-  function(x, family, alpha, theta){
  if(ncol(as.matrix(x))>1) stop("x should be a vector or a n by 1 matrix")
  alpha <- as.integer(alpha)
  if(alpha==0){stop("Invalid parameter alpha")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  val <-  switch(family, 
                gumbel = x^(alpha-1)/gamma(alpha)*copGumbel@absdPsi(x, theta = theta, degree = alpha), 
                clayton = x^(alpha-1)/gamma(alpha)*copClayton@absdPsi(t = theta*x, theta, degree = alpha)*(theta^alpha), 
                frank = x^(alpha-1)/gamma(alpha)*copFrank@absdPsi(x, theta = theta, degree = alpha), 
                AMH = x^(alpha-1)/gamma(alpha)*copAMH@absdPsi(x, theta = theta, degree = alpha), 
                joe = x^(alpha-1)/gamma(alpha)*copJoe@absdPsi(x, theta = theta, degree = alpha))
  return(pmax(val, 0))
}


##########################################################################
#' Old code for the copula likelihood function for Liouville copulas
#'
#' The function is used internally for optimization.
#'
#' @param theta parameter of the corresponding Archimedean copula
#' @param data sample matrix from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param MC.approx whether to use Monte-Carlo approximation for the inverse survival function (default is \code{TRUE})
#'
#' @return value of marginal density
.pliouv.opt_old <-  function(theta, data, family, alphavec, MC.approx = T){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  #changed order of parameters for optimization functions
  illegalpar <-  switch(family, 
                       clayton = copClayton@paraConstr(theta), 
                       gumbel = copGumbel@paraConstr(theta), 
                       frank = copFrank@paraConstr(theta), 
                       AMH = copAMH@paraConstr(theta), 
                       joe = copJoe@paraConstr(theta))
  if(!illegalpar)
    stop("Illegal parameter value")
  n <-  dim(data)[1]
  d <-  dim(data)[2]
  Xdata = data
  if(MC.approx == T){
    Xdata <-  H_inv(data, alphavec = alphavec, family = family, theta = theta)
  }
  Xmargin <-  data
  for (j in 1:d)
  {
    if(MC.approx == F){
      Xdata[, j] <-  isliouvm(data[, j], family, alphavec[j], theta)
    }
    Xmargin[, j] <-  log(dliouvm(Xdata[, j], family, alphavec[j], theta))
  }
  numerator.log <-  dliouv(Xdata, family, alphavec, theta, is.log = T)
  denominator.log <-  sum(Xmargin)
  numerator.log - denominator.log
}

##########################################################################
#' Internal code for the copula likelihood function for Liouville copulas
#'
#' The function is used internally for optimization.
#'
#' @param theta parameter of the corresponding Archimedean copula
#' @param data sample matrix from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param MC.approx whether to use Monte-Carlo approximation for the inverse survival function (default is \code{TRUE})
#'
#' @return value of marginal density
pliouv.opt <-  function(theta, data, family, alphavec, MC.approx = T){
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  #changed order of parameters for optimization functions
  illegalpar <-  switch(family, 
                       clayton = copClayton@paraConstr(theta), 
                       gumbel = copGumbel@paraConstr(theta), 
                       frank = copFrank@paraConstr(theta), 
                       AMH = copAMH@paraConstr(theta), 
                       joe = copJoe@paraConstr(theta))
  if(!illegalpar)
    stop("Illegal parameter value")
  n <-  dim(data)[1]
  d <-  dim(data)[2]
  Xdata = data
  if(MC.approx == T){
    Xdata <-  H_inv(data, alphavec = alphavec, family = family, theta = theta)
  }
  Xmargin <-  data
  for (j in 1:d){
    if(MC.approx == F){
      Xdata[, j] <-  isliouvm(data[, j], family, alphavec[j], theta)
    }
    Xmargin[, j] <-  switch(family, 
                          gumbel = copGumbel@absdPsi(Xdata[, j], theta = theta, degree = alphavec[j], log = T), 
                          clayton = copClayton@absdPsi(t = theta*Xdata[, j], theta, degree = alphavec[j], log = T), 
                          frank = copFrank@absdPsi(Xdata[, j], theta = theta, degree = alphavec[j], log = T), 
                          AMH = copAMH@absdPsi(Xdata[, j], theta = theta, degree = alphavec[j], log = T), 
                          joe = copJoe@absdPsi(Xdata[, j], theta = theta, degree = alphavec[j], log = T)
    )
  }
  alpha = sum(alphavec)
  numerator.log <-  sum(switch(family, 
                              gumbel = copGumbel@absdPsi(rowSums(Xdata), theta = theta, degree = alpha, log = T), 
                              clayton = copClayton@absdPsi(t = theta*rowSums(Xdata), theta, degree = alpha, log = T), 
                              frank = copFrank@absdPsi(rowSums(Xdata), theta = theta, degree = alpha, log = T), 
                              AMH = copAMH@absdPsi(rowSums(Xdata), theta = theta, degree = alpha, log = T), 
                              joe = copJoe@absdPsi(rowSums(Xdata), theta = theta, degree = alpha, log = T)
  ))
  denominator.log <-  sum(Xmargin)
  numerator.log - denominator.log
}
##########################################################################
#' Maximize copula likelihood function for Liouville copulas
#'
#' A wrapper to \code{optim} using the Nelder-Mead algorithm to maxime pointwise given every alphavec over a grid. 
#' Returns the maximum for \code{alphavec} and \code{theta}.
#'
#' @param data sample matrix from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param interval interval over which to look for \code{theta} (bounds for Nelder-Mead)
#' @param boundary vector of endpoints for search of Dirichlet allocation parameters. Either \code{boundary} or \code{lattice.mat} can be supplied
#' @param lattice.mat matrix of tuples of Dirichlet allocation parameters at which to evaluate the likelihood
#' @param return_all should all results (as list) or only maximum value be returned. Defaults to \code{FALSE}
#' @param MC.approx whether to use Monte-Carlo approximation for the inverse survival function (default is \code{TRUE})
#'
#' @return a list with values of \code{theta} and Dirichlet parameter along with maximum found. Gives index of maximum amongst models fitted.
#' @export
#' @examples 
#' ## data <- rliouv(n=100, family="joe", alphavec=c(1,2), theta=2)
#' ## liouv.maxim(data=data, family="j", interval=c(1.25,3), boundary=c(2,2),return_all=TRUE)
#' ## lattice.mat <- t(combn(1:3,2))
#' ## liouv.maxim(data=data, family="j", interval=c(1.25,3), lattice.mat=lattice.mat, return_all=FALSE)
liouv.maxim <- function(data, family, interval, boundary = NULL, lattice.mat = NULL, return_all = FALSE, MC.approx = TRUE){
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  stopifnot(!is.null(lattice.mat) || !is.null(boundary))
  if(!is.null(boundary)){stopifnot(length(boundary) == dim(data)[2]) }
  if(!is.null(lattice.mat)){stopifnot(ncol(lattice.mat) == ncol(data)) }
    stopifnot(length(interval) == 2, interval[1]<interval[2], 
            switch(family, clayton = copClayton@paraConstr(interval[1]), 
                     gumbel = copGumbel@paraConstr(interval[1]), 
                     frank = copFrank@paraConstr(interval[1]), 
                     AMH = copAMH@paraConstr(interval[1]), 
                     joe = copJoe@paraConstr(interval[1])), 
              switch(family, clayton = copClayton@paraConstr(interval[2]), 
                     gumbel = copGumbel@paraConstr(interval[2]), 
                     frank = copFrank@paraConstr(interval[2]), 
                     AMH = copAMH@paraConstr(interval[2]), 
                     joe = copJoe@paraConstr(interval[2])))
  
  if(is.null(lattice.mat)){
    term <- list()
    for(i in 1:length(boundary)){
      term[[i]] = 1:boundary[i]
    }
    lattice.mat <- matrix(unlist(expand.grid(term)),ncol=length(boundary))
  }
  #Define containers
  like_array  <- array(dim=apply(lattice.mat,2,max))
  theta_array <- array(dim=apply(lattice.mat,2,max))
  #Optimization over theta
  for(i in 1:nrow(lattice.mat)){
  lattice.vec <- lattice.mat[i,]
  # Choice of maximization function, see http://cran.r-project.org/web/views/Optimization.html
  res <- optimise(pliouv.opt, data = data, family = family, alphavec = lattice.vec, 
             interval = interval, maximum = T, MC.approx = MC.approx)
   theta_array[matrix(lattice.vec, 1)] <- res$maximum
    like_array[matrix(lattice.vec, 1)] <- res$objective
  }
  indexed <- which.max(like_array)
  pos <- which( like_array==max(like_array,na.rm=T) , arr.ind = T )
  if(return_all){
    return(list(index = pos, 
                theta_max = theta_array[indexed], 
                theta_val = theta_array, 
                log.lik_max = like_array[indexed], 
                log.lik = like_array))
  }  else{
    return(list(index = pos, 
                theta_max = theta_array[indexed], log.lik = like_array[indexed]))
  }
}
##########################################################################
#' Maximize copula likelihood function for Liouville copulas via methods of moments
#'
#' A wrapper to \code{optim} using the methods of moments to maxime pointwise given every \code{alphavec} over a grid. 
#' Returns the maximum for \code{alphavec} and \code{theta}
#'
#' @param data sample matrix from a Liouville copula
#' @param family family of the Liouville copula. Either \code{"clayton"} or \code{"gumbel"}
#' @param boundary vector of endpoints for search of Dirichlet allocation parameters. Either \code{boundary} or \code{lattice.mat} can be supplied
#' @param lattice.mat matrix of tuples of Dirichlet allocation parameters at which to evaluate the likelihood
#' @param return_all should all results (as list) or only maximum value be returned. Defaults to \code{FALSE}
#' @param MC.approx whether to use Monte-Carlo approximation for the inverse survival function (default is \code{TRUE})
#'
#' @return a list with values of theta and Dirichlet and maximum found. Gives index
#' @export
#' @examples
#' data <- rliouv(n=1000, family="gumbel", alphavec=c(1,2), theta=2)
#' liouv.maxim.mm(data=data, family="gumbel", boundary=c(3,3),return_all=TRUE)
#' lattice.mat <- t(combn(1:3,2))
#' liouv.maxim.mm(data=data, family="gumbel", lattice.mat=lattice.mat, return_all=FALSE)
liouv.maxim.mm <- function(data, family, boundary = NULL, lattice.mat = NULL, return_all = FALSE, MC.approx = T){
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  stopifnot(!is.null(lattice.mat) || !is.null(boundary))
  if(is.null(lattice.mat)){
    stopifnot(length(boundary) == dim(data)[2] && dim(data)[2] == 2)
    term <- list()
    for(i in 1:length(boundary)){
      term[[i]] = 1:as.integer(boundary[i])
    }
    lattice.mat <- matrix(unlist(expand.grid(term)),ncol=length(boundary))
  }
  if(!is.null(lattice.mat)){stopifnot(ncol(lattice.mat) == ncol(data), ncol(data) == 2) }
  tau_hat <- cor.fk(data)[1, 2]
  #Create containers for values
  like_array  <- array(dim=apply(lattice.mat,2,max))
  theta_array <- array(dim=apply(lattice.mat,2,max))
  #Optimization over theta
  for(i in 1:nrow(lattice.mat)){
  lattice.vec <- lattice.mat[i,]
    theta_hat<- liouv.iTau(tau_hat, family, alphavec = lattice.vec)
    theta_array[matrix(lattice.vec, 1)] <- theta_hat
    like_array[matrix(lattice.vec, 1)] <- pliouv.opt(theta_hat, data, family, lattice.vec, MC.approx = MC.approx)
   
  }
  #Which is the max
  indexed <- which.max(like_array)
  pos <- which( like_array==max(like_array,na.rm=T) , arr.ind = T )
  if(return_all){
    return(list(index = pos, 
                theta_max = theta_array[indexed], 
                theta_val = theta_array, 
                log.lik_max = like_array[indexed], 
                log.lik = like_array))
  }  else{
    return(list(index = pos, 
                theta_max = theta_array[indexed], log.lik = like_array[indexed]))
  }
}

##########################################################################
#' Computes Kendall's \eqn{\tau}{tau} for Clayton or Gumbel Liouville copula
#'
#' The function computes Kendall's \eqn{\tau}{tau} for the given model, given \code{alphavec}
#'
#' @param theta parameter of the corresponding Archimedean copula
#' @param family family of the Liouville copula. Either \code{"clayton"} or \code{"gumbel"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#'
#' @return value of \eqn{\tau}{tau}
.liouv.Tau_s <- function(theta, family, alphavec){
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  if(length(alphavec) != 2)
    stop("Illegal parameter length - only implemented for the bivariate case")
  a1 <- alphavec[1]; a2 <- alphavec[2]; a <- a1+a2

  if(family == "clayton"){
    express_coef_clayton <- function(i, j, theta, a){
      return(exp(lgamma(1/theta+i+j)+lgamma(1/theta+a)+lgamma(2/theta)+lgamma(a+i+j)-
                   2*log(theta)-2*lgamma(1/theta+1)-lgamma(2/theta+i+j+a)-lgamma(a)))
    }
  }
  if(family == "gumbel"){
    express_coef_gumb <- function(index_k, index_l, theta, a){
      if(index_l == 0){
        d1 <- .coeffG(d = index_k, alpha = 1/theta, method = "direct", log = T)
        intern <- seq(1:index_k)
        d2 <- lgamma(intern)-intern*log(2)
        return(theta*sum(exp(d1+d2))/gamma(a))
      }      else{
        d1 <- outer(.coeffG(d = index_k, alpha = 1/theta, method = "direct", log = T), 
                  .coeffG(d = index_l, alpha = 1/theta, method = "direct", log = T), "+")
        intern <- outer(seq(1:index_k), seq(1:index_l), "+")
        d2 <- lgamma(intern)-intern*log(2)
        return(theta*sum(exp(d1+d2))/gamma(a))
      }
    }
  }
  count <- sum(sapply(0:(a1-1), function(i){
    sum(sapply(0:(a2-1), function(j){
      exp(lbeta(a1+i, a2+j)-lbeta(a1, a2)-lgamma(i+1)-lgamma(j+1))*
        switch(family, 
               gumbel = express_coef_gumb(index_k = a, index_l = i+j, theta, a = a), 
               clayton = express_coef_clayton(i = i, j = j, theta, a = a))
    }))
  }))
  return(count*4-1)
}
#' Computes Kendall's tau for Clayton or Gumbel Liouville copula
#'
#' The function computes Kendall's \eqn{\tau}{tau} for the given model, given \code{alphavec}
#'
#' @param theta parameter of the corresponding Archimedean copula
#' @param family family of the Liouville copula. Either \code{"clayton"} or \code{"gumbel"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @examples
#' liouv.Tau(theta=2, family="gumbel", alphavec=c(1,2))
#' liouv.Tau(theta=1, family="clayton", alphavec=c(2,1))
#' @return vector of \eqn{\tau}{tau}
#' @export
liouv.Tau = Vectorize(FUN = .liouv.Tau_s, vectorize.args = "theta")

##########################################################################
#' Moment estimate for \code{theta} derived from Kendall's \eqn{\tau}{tau}
#'
#' The moment estimates are based on inversion of the formula Kendall's \eqn{\tau}{tau} for Clayton or Gumbel Liouville copula
#'
#' @param tau_hat estimated Kendall's \eqn{\tau}{tau} value from data
#' @param family family of the Liouville copula. Either \code{"clayton"} or \code{"gumbel"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#'
#' @return Value of theta
.liouv.iTau_s <-  function(tau_hat, family, alphavec){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel"))
  if(family != "gumbel" && family != "clayton")
    stop("Not yet implemented for this family")
  invtau <-  function(theta, tau, family, alphavec)
  {
    liouv.Tau(theta, family, alphavec)-tau
  }
  lb <- switch(family, gumbel = 1+1e-60, clayton = 1e-60)
  #Not implemented for negative theta values (Clayton family)
  out <-  uniroot(invtau, lower = lb, upper = 100, f.lower = -1, f.upper = 1, 
                 alphavec = alphavec, tau = tau_hat, family = family)
  out$root
}
#' Moment estimate for theta derived from Kendall's tau formula
#'
#' The moment estimates are based on inversion of the formula Kendall's tau for Clayton or Gumbel Liouville copula
#'
#' @param tau_hat estimated vector of Kendall's tau values
#' @param family family of the Liouville copula. Either \code{"clayton"} or \code{"gumbel"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @examples 
#' liouv.iTau(0.5,family="gumbel", c(1,2))
#' liouv.iTau(0.5,family="clayton", c(3,2))
#' @return Vector of theta
#' @export
liouv.iTau = Vectorize(.liouv.iTau_s, "tau_hat")

##########################################################################
#' Inverse survival function if Monte-Carlo approximation is set to \code{TRUE} in \code{liouv.maxim}
#'
#' The function is used internally for optimization.
#'
#' @keywords internal
#'
#' @param u data at which to compute the survival inverse
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param theta parameter of the corresponding Archimedean copula
#' @param MC number of Monte-Carlo points for evaluation
#' @param trunc whether to truncate at low quantile. Will be based on numerical root finding for the lower 0.025 fraction of the data
#'
#' @return Inverse survival function values
#' @examples
#' u <- rliouv(n = 100, family = "frank", alphavec <- c(2,3), theta = 1)
#' H_inv(u=u, family="frank", alphavec=c(2,3), theta=2)
#' #Difference between true value and approximation (can be large depending on family)
#' sum(abs(H_inv(u=u, family="frank", alphavec=c(2,3), theta=2)-
#' isliouvm_m(u=u, family="frank", alphavec=c(2,3), theta=2)))
H_inv <-  function(u, alphavec, family, theta, MC = 100000, TRUNC = F){
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  d <-  sum(alphavec)
  #Parametrization for Clayton is different
  if(family == "clayton"){
    mc <- 	(rCopula(n = MC, claytonCopula(param = theta, dim = d))^(-theta)-1)/theta
  }	else if(family=="AMH") {
     mc <-  iPsi(get(paste(tolower(family), "Copula", sep = ""))(param = theta), 
               rarchi(n = MC, family="AMH", theta = theta, d = d))
  } else{
    mc <-  iPsi(get(paste(tolower(family), "Copula", sep = ""))(param = theta), 
               rCopula(n = MC, get(paste(tolower(family), "Copula", sep = ""))(param = theta, dim = d)))
  }
  mc <-  .marginCombo(alphavec, mc)
  inverse = numeric()
  for(i in 1:length(alphavec)){
    inverse = cbind(inverse, as.numeric(quantile(mc[, i], probs = (1-u[, i]))))
  }
  if(TRUNC == T){
    hquant = which(u<0.025, arr.ind = T)
    for(i in 1:(dim(hquant)[1])){
      inverse[hquant[i, 1], hquant[i, 2]] = 	isliouvm(u[hquant[i, 1], hquant[i, 2]], family, alphavec[hquant[i, 2]], theta)
    }

    hquant = which(u>0.975, arr.ind = T)
    for(i in 1:(dim(hquant)[1])){
      inverse[hquant[i, 1], hquant[i, 2]] = 	isliouvm(u[hquant[i, 1], hquant[i, 2]], family, alphavec[hquant[i, 2]], theta)
    }
  }
  return(inverse)
}
##########################################################################
### Rcpp function to sum margins - i.e. rowSums
#marginCombo() in Rcpp file
##########################################################################
#' Clayton coefficient for Kendall's tau formula
#'
#' The function is used internally in \cite{liouv.Tau}.
#'
#' @keywords internal
#'
#' @param i first index
#' @param j second index
#' @param theta parameter of the corresponding Archimedean copula
#' @param a \eqn{L_1}{L1}-norm of Dirichlet allocation vector
#'
#' @return Clayton coefficient
express_coef_clayton <- function(i, j, theta, a){
  return(exp(lgamma(1/theta+i+j)+lgamma(1/theta+a)+lgamma(2/theta)+lgamma(a+i+j)-
               2*log(theta)-2*lgamma(1/theta+1)-lgamma(2/theta+i+j+a)-lgamma(a)))
}

# Note: function coeffG is a copy of that found in the `copula' package
# the latter is unexported, and CRAN won't allow access via :::
# All copyrights and credits to copula authors
# (c) Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
.coeffG <-  function (d, alpha, method = c("sort", "horner", "direct", "dsumSibuya", 
paste("dsSib", eval(formals(dsumSibuya)$method), sep = ".")), 
log = FALSE, verbose = FALSE) 
{
stopifnot(is.numeric(d), length(d) == 1, d >= 1, length(alpha) == 
    1, 0 < alpha, alpha <= 1)
a <-  numeric(d)
method <-  match.arg(method)
if (method == "dsumSibuya") {
    message("method 'dsumSibuya' is deprecated.\n  ", "use method = 'dsSib.log' instead")
    method <-  "dsSib.log"
}
Meth <-  if (grepl("^dsSib", method)) {
    meth.dsSib <-  sub("^dsSib\\.(.*)", "\\1", method)
    "dsSib"
}
else method
switch(Meth, sort = {
    ls <-  log(abs(Stirling1.all(d)))
    lS <-  lapply(1:d, function(n) log(Stirling2.all(n)))
    wrong.sign <-  integer()
    for (k in 1:d) {
        j <-  k:d
        b <-  j * log(alpha) + ls[j] + unlist(lapply(j, function(i) lS[[i]][k]))
        b.max <-  max(b)
        exponents <-  b - b.max
        exps <-  exp(exponents)
        even <-  if (k == d) NULL else seq(2, d - k + 1, by = 2)
        odd <-  seq(1, d - k + 1, by = 2)
        sum.neg <-  sum(sort(exps[even]))
        sum.pos <-  sum(sort(exps[odd]))
        sum. <-  sum.pos - sum.neg
        a[k] <-  if (log) b.max + log(sum.) else exp(b.max) * 
            sum.
        if (sum.neg > sum.pos) {
            if (verbose) message("sum.neg > sum.pos for k = ", 
              k)
            wrong.sign <-  c(wrong.sign, k)
        }
    }
    if (length(wrong.sign)) attr(a, "wrong.signs") <-  wrong.sign
    a
}, dsSib = {
    k <-  1:d
    ck <-  if (log) c(0, cumsum(log(d:2)))[d:1] else c(1, 
        cumprod(d:2))[d:1]
    p <-  dsumSibuya(d, k, alpha, method = meth.dsSib, log = log)
    if (log) p + ck else p * ck
}, horner = {
    s.abs <-  abs(Stirling1.all(d))
    k <-  1:d
    S <-  lapply(k, Stirling2.all)
    pol <-  vapply(k, function(k.) {
        j <-  0:(d - k.)
        c.j <-  s.abs[j + k.] * vapply(j, function(i) S[[i + 
            k.]][k.], 1)
        polynEval(c.j, -alpha)
    }, NA_real_)
    if (log) k * log(alpha) + log(pol) else alpha^k * pol
}, direct = {
    s <-  Stirling1.all(d)
    k <-  1:d
    S <-  lapply(k, Stirling2.all)
    vapply(k, function(k.) {
        j <-  k.:d
        S. <-  vapply(j, function(i) S[[i]][k.], 1)
        sm <-  sum(alpha^j * s[j] * S.)
        if (log) log(abs(sm)) else (-1)^(d - k.) * sm
    }, NA_real_)
}, stop(gettextf("unsupported method '%s' in coeffG", method), 
    domain = NA))
}

##########################################################################
#' Gumbel coefficient for Kendall's tau formula
#'
#' The function is used internally in \cite{liouv.Tau}.
#'
#' @keywords internal
#'
#' @param index_k first index
#' @param index_l second index
#' @param theta parameter of the corresponding Archimedean copula
#' @param a sum of alphavec
#'
#' @return Gumbel coefficient
express_coef_gumb <- function(index_k, index_l, theta, a){
  if(index_l == 0){
    d1 <- .coeffG(d = index_k, alpha = 1/theta, method = "direct", log = T)
    intern <- seq(1:index_k)
    d2 <- lgamma(intern)-intern*log(2)
    return(theta*sum(exp(d1+d2))/gamma(a))
  }  else{
    d1 <- outer(.coeffG(d = index_k, alpha = 1/theta, method = "direct", log = T), 
              .coeffG(d = index_l, alpha = 1/theta, method = "direct", log = T), "+")
    intern <- outer(seq(1:index_k), seq(1:index_l), "+")
    d2 <- lgamma(intern)-intern*log(2)
    return(theta*sum(exp(d1+d2))/gamma(a))
  }
}


##########################################################################
#' Parametric bootstrap confidence interval for the parameter \code{theta} for Liouville copula
#'
#' The parametric bootstrap provides confidence intervals by repeatedly sampling datasets from the postulated 
#' Liouvilla copula model. If \eqn{d=2} and the model is either \code{gumbel} or \code{clayton}, the value of
#' Kendall's \eqn{\tau}{tau} is calculated from the sample, and the confidence interval or the quantiles correspond
#' to the inverse \eqn{\tau^{-1}(\tau(\theta))}{tau} for the bootstrap quantile values of \eqn{\tau}{tau} (using monotonicity).
#' 
#' Since no closed-form formulas exist for the other models or in higher dimension, 
#' the method is extremely slow since it relies on maximization 
#' of a new sample from the model and look up the corresponding parameters.
#'
#' @param B number of bootstrap replicates
#' @param family family of the Liouville copula. Either \code{"clayton"}, \code{"gumbel"}, \code{"frank"}, \code{"AMH"} or \code{"joe"}
#' @param alphavec vector of Dirichlet allocations (must be a vector of integers)
#' @param n sample size
#' @param theta.hat estimate of theta
#' @param quant if the vector of probability is specified, the function will return the corresponding bootstrap quantiles
#' @param silent boolean for output progress. Default is \code{FALSE}, which means iterations are printed if \eqn{d>2}.
#' @examples 
#' theta.bci(B=99, family="gumbel", alphavec=c(2,3), n=1000, theta.hat=2)
#' ## theta.bci(B=19, family="AMH", alphavec=c(1,2), n=100, theta.hat=0.5, quant=c(0.05,0.95))
#' ## theta.bci(B=19, family="frank", alphavec=c(1,2,3), n=100, theta.hat=0.5, quant=c(0.05,0.95))
#' @return a list with a 95% confidence inteval unless selected quantiles \code{quant} are supplied 
#' and the bootstrap values of Kendall's tau in \code{boot_tau} if \eqn{d=2} and the model is either \code{gumbel} or \code{clayton}.
#' Otherwise, the list contains \code{boot_theta}.
#' @export
theta.bci <- function(B = 1999, family, alphavec, n, theta.hat, quant = c(0.025,0.975), silent=F){
  family <-  match.arg(family, c("clayton", "gumbel", "frank", "AMH", "joe"))
  alphavec <- as.integer(alphavec)
  if(any(alphavec)==0){stop("Invalid parameters")}
  illegalpar <-  switch(family, 
                       clayton = copClayton@paraConstr(theta.hat), 
                       gumbel = copGumbel@paraConstr(theta.hat), 
                       frank = copFrank@paraConstr(theta.hat), 
                       AMH = copAMH@paraConstr(theta.hat), 
                       joe = copJoe@paraConstr(theta.hat))
  if(!illegalpar)
    stop("Illegal parameter value")
  d <- length(alphavec)
  if(d == 2 && family %in% c("gumbel","clayton")){
    boot.rep <- sapply(1:B, function(x){
    cor.fk(rliouv(n = n, family, alphavec, theta.hat))[1, 2]
    })
    quantiles = NaN
    if(is.vector(quant)){
      quantiles = liouv.iTau(quantile(boot.rep, probs = quant), family, alphavec)
    }
    return(list(quantiles = quantiles, boot_kendall = boot.rep))
  }  else{
    boot.rep <- sapply(1:B, function(x){
      if(silent==FALSE){print(paste("Step", x))}
      boot.samp <- rliouv(n = n, family, alphavec, theta.hat)
      boot.opt <- liouv.maxim(boot.samp, family, lattice.mat = as.matrix(t(alphavec)), MC.approx = F, 
                            interval = switch(family, gumbel = c(1.001, 10), 
                                            clayton = c(0.001, 10), 
                                            AMH = c(0.001, 0.999), 
                                            joe = c(1.001, 10), 
                                            frank = c(0.001, 10)))
      boot.opt$theta_max})
    quantiles = NaN
    if(is.vector(quant)){
      quantiles = quantile(boot.rep, probs = quant)
    }
    return(list(quantiles = quantiles, boot_theta = boot.rep))

  }
}
