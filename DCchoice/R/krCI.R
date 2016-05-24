krCI <- function(obj = NULL, nsim = 1000, CI = 0.95, individual = NULL){
  if(CI > 1) stop("CI must be between 0 and 1")
  if(class(obj) != "sbchoice" & class(obj) != "dbchoice"){
    # stop if the object is neither a sdchoice nor a dbchoice class
    stop("the object must be either dbchoice of sbchoice class")
  }
  
  X <- obj$covariates   # retrieving the covariates from the object
  npar <- length(obj$coefficients)   # the number of coefficients
  coef <- obj$coefficients   # coefficient estimates

  formula <- formula(obj$formula, lhs = 0, rhs = -2)

  if(class(obj) == "sbchoice"){   # covariance matrix of the estimates
    if(obj$dist != "weibull"){
      V <- vcov(obj$glm.out)
    } else {
      V <- solve(obj$glm.out$hessian)
    }
  } else {
    V <- solve(obj$Hessian)
  }

  kr.coef <- t(mvrnorm(nsim, coef, V))    # mvnorm is from MASS package

  kr.b <- kr.coef[npar, ]                   # coefficients on the covariates
  kr.coef <- kr.coef[-npar, , drop = FALSE]               # coefficient on the bid
  if(is.null(individual)) {
    kr.Xb <- colSums(colMeans(X[, -npar, drop = FALSE])*kr.coef)
  } else {
    mf.nexX <- model.frame(formula, individual, xlev = obj$xlevels)
    mm.newX <- model.matrix(formula, mf.nexX, contrasts.arg = obj$contrasts)
    newX <- as.vector(mm.newX)
    kr.Xb <- colSums(newX * kr.coef)
  }

  if(obj$dist == "log-logistic"){ # for log-logistic error distribution
    mbid <- exp(max(obj$bid))     # the maximum bid
    kr.median <- sort(exp(-kr.Xb/kr.b))
    func.kr <- function(x, A, B) plogis(-(A + B*log(x)), lower.tail = FALSE)
    denom <- function(A, B) plogis(-(A + B*log(mbid)))
  } else if(obj$dist == "log-normal"){  # for log-normal error distribution
    mbid <- exp(max(obj$bid))
    kr.median <- sort(exp(-kr.Xb/kr.b))
    func.kr <- function(x, A, B) pnorm(-(A + B*log(x)), lower.tail = FALSE)
    denom <- function(A, B) pnorm(-(A + B*log(mbid)))
  } else if(obj$dist == "logistic") {  # for logistic error distribution
    mbid <- max(obj$bid)
    kr.median <- sort(-kr.Xb/kr.b)
    func.kr <- function(x, A, B) plogis(-(A + B*x), lower.tail = FALSE)
    denom <- function(A, B) plogis(-(A + B*mbid))
  } else if(obj$dist == "normal") {  # for normal error distribution
    mbid <- max(obj$bid)
    kr.median <- sort(-kr.Xb/kr.b)
    func.kr <- function(x, A, B) pnorm(-(A + B*x), lower.tail = FALSE)
    denom <- function(A, B) pnorm(-(A + B*mbid))
  } else if(obj$dist == "weibull") {
    mbid <- exp(max(obj$bid))
    kr.median <- sort(exp(-kr.Xb/kr.b)*(log(2))^(-1/kr.b))
    func.kr <- function(x, A, B) pweibull(exp(-A - B*log(x)), shape=1, lower.tail=FALSE)
    denom <- function(A, B) pweibull(exp(-A - B*log(mbid)), shape=1)
  }

  kr.meanWTP <- numeric(nsim)
  kr.trunc.meanWTP <- numeric(nsim)
  kr.adj.trunc.meanWTP <- numeric(nsim)
  
  for(i in 1:nsim){
    kr.meanWTP[i] <- integrate(func.kr, 0, Inf, A = kr.Xb[i], B = kr.b[i], stop.on.error = FALSE)$value
    kr.trunc.meanWTP[i] <- integrate(func.kr, 0, mbid, A = kr.Xb[i], B = kr.b[i], stop.on.error = FALSE)$value
    kr.adj.trunc.meanWTP[i] <- integrate(func.kr, 0, mbid, A = kr.Xb[i], B = kr.b[i], stop.on.error = FALSE)$value/denom(A = kr.Xb[i], B = kr.b[i])
  }
  
  output <- list(mWTP = kr.meanWTP, tr.mWTP = kr.trunc.meanWTP, adj.tr.mWTP = kr.adj.trunc.meanWTP, medWTP = kr.median)
  
  kr.meanWTP <- sort(kr.meanWTP)
  kr.trunc.meanWTP <- sort(kr.trunc.meanWTP)
  kr.adj.trunc.meanWTP <- sort(kr.adj.trunc.meanWTP)

 lb <- 0.5*(1 - CI)   # lower bound of the empirical distribution
 ub <- CI + lb        # upper bound of the empirical distribution

  int <- c(ceiling(nsim*lb), floor(nsim*ub))  # subscripts corresponding to lb and ub

  # 100*CI% simulated conficdence intervals
  CImat <- rbind(kr.meanWTP[int],       # for mean
                 kr.trunc.meanWTP[int], # for truncated mean
                 kr.adj.trunc.meanWTP[int], # for truncated mean with adjustment
                 kr.median[int])        # for median

  if (is.null(individual)) {
    tmp.sum <- wtp(object = obj$covariates, b = obj$coefficients, bid = obj$bid, dist = obj$dist)
  } else {
    tmp.sum <- wtp(object = mm.newX, b = obj$coefficients, bid = obj$bid, dist = obj$dist)
  }

#  sim.mean <- c(mean(kr.meanWTP), mean(kr.trunc.meanWTP), mean(kr.adj.trunc.meanWTP), median(kr.median))
#  sim.se <- c(sd(kr.meanWTP), sd(kr.trunc.meanWTP), sd(kr.adj.trunc.meanWTP), -999)
  
  out <- cbind(c(tmp.sum$meanWTP, tmp.sum$trunc.meanWTP, tmp.sum$adj.trunc.meanWTP, tmp.sum$medianWTP), 
               CImat)
  rownames(out) <- c("Mean", "truncated Mean", "adjusted truncated Mean", "Median")
  colnames(out) <- c("Estimate", "LB", "UB")

  if(!is.finite(tmp.sum$meanWTP)){  # when the parameter estimate does not satisfy the finite mean WTP condition
    out[1, 2:3] <- -999     # set the interval [-999, -999]
    output$mWTP <- -999      # set the simulated mean to NULL
  }  

  output$out <- out
  
  class(output) <- "krCI"
  return(output)
  
}

print.krCI <- function(x, ...){
  cat("the Krinsky and Robb simulated confidence intervals\n")
  printCoefmat(x$out, digits = 5)
  invisible(x$out)
}

# summary.krCI <- function(object, ...){
#   
#   invisible(object)
# }


