# Computing a bootsrap confidence interval
bootCI <- function (obj, nboot = 1000, CI = 0.95, individual = NULL){
  if(class(obj) != "sbchoice" & class(obj) != "dbchoice"){
    # stop if the object is neither a sdchoice nor a dbchoice class
    stop("the object must be either dbchoice of sbchoice class")
  }
  
  if(CI > 1 | CI < 0) stop("CI must be between 0 and 1")
  
    tmp.dat <- eval(obj$data.name, parent.frame())   # retrieving the data from the object

  nobs <- obj$nobs
  ind <- 1:nobs
  boot.mean <- numeric(nboot)
  boot.median <- numeric(nboot)
  boot.trmean <- numeric(nboot)
  boot.adj.trmean <- numeric(nboot)

  dist <- obj$distribution
  fr <- obj$formula

  if (!is.null(individual)) {
    formula <- formula(obj$formula, lhs = 0, rhs = -2)
    mf.nexX <- model.frame(formula, individual, xlev = obj$xlevels)
    mm.newX <- model.matrix(formula, mf.nexX, contrasts.arg = obj$contrasts)
  }

  if(class(obj) == "dbchoice"){
    for(i in 1:nboot){
      ind.boot <- sample(ind, nobs, replace = TRUE)  # determining the number of rows for bootstrap sample with replacement
      boot.dat <- tmp.dat[ind.boot, ]  # resampling data
      suppressWarnings(
        tmp.re <- dbchoice(fr, data = boot.dat, dist = dist, par=obj$coefficients)  # estimating a DBCV model
      )
      if(tmp.re$convergence){
        if (is.null(individual)) {
          boot <- wtp(object = tmp.re$covariates, b = tmp.re$coefficients, bid = tmp.re$bid, dist = tmp.re$dist)
        } else {
          boot <- wtp(object = mm.newX, b = tmp.re$coefficients, bid = tmp.re$bid, dist = tmp.re$dist)
        }
        boot.mean[i] <- boot$meanWTP
        boot.median[i] <- boot$medianWTP
        boot.trmean[i] <- boot$trunc.meanWTP
        boot.adj.trmean[i] <- boot$adj.trunc.meanWTP
      } else {
        i < i - 1        # discard an unconverged trial
      }
    }
  } else if(class(obj) == "sbchoice"){
    for(i in 1:nboot){
      ind.boot <- sample(ind, nobs, replace = TRUE)
      boot.dat <- tmp.dat[ind.boot, ]
      suppressWarnings(
        tmp.re <- sbchoice(fr, data = boot.dat, dist = dist, par=obj$coefficients)
      )
      if(tmp.re$glm.out$converged){
        if (is.null(individual)) {
          boot <- wtp(object = tmp.re$covariates, b = tmp.re$coefficients, bid = tmp.re$bid, dist = tmp.re$dist)
        } else {
          boot <- wtp(object = mm.newX, b = tmp.re$coefficients, bid = tmp.re$bid, dist = tmp.re$dist)
        }
        boot.mean[i] <- boot$meanWTP
        boot.median[i] <- boot$medianWTP
        boot.trmean[i] <- boot$trunc.meanWTP
        boot.adj.trmean[i] <- boot$adj.trunc.meanWTP
      } else {
        i < i -1        # discard an unconverged trial
      }
    }
  }
  
  output <- list(mWTP = boot.mean, tr.mWTP = boot.trmean, 
                adj.tr.mWTP = boot.adj.trmean, medWTP = boot.median)
  
  # sorting the simulation outcomes
  boot.mean <- sort(boot.mean)
  boot.median <- sort(boot.median)
  boot.trmean <- sort(boot.trmean)
  boot.adj.trmean <- sort(boot.adj.trmean)

  lb <- 0.5*(1 - CI)   # lower bound of the empirical distribution
  ub <- CI + lb        # upper bound of the empirical distribution

  int <- c(ceiling(nboot*lb), floor(nboot*ub))  # subscripts corresponding to lb and ub

  # 100*CI% bootstrap CI
  CImat <- rbind(boot.mean[int],    # for mean
                 boot.trmean[int],  # for truncated mean
                 boot.adj.trmean[int],  # for truncated mean with adjustment
                 boot.median[int])  # for median
  
  # the mean estimates in the original outcome
  if (is.null(individual)) {
    tmp.sum <- wtp(object = obj$covariates, b = obj$coefficients, bid = obj$bid, dist = obj$dist)
  } else {
    tmp.sum <- wtp(object = mm.newX, b = obj$coefficients, bid = obj$bid, dist = obj$dist)
  }
  
  out <- cbind(c(tmp.sum$meanWTP, tmp.sum$trunc.meanWTP, tmp.sum$adj.trunc.meanWTP, tmp.sum$medianWTP), CImat)
  rownames(out) <- c("Mean", "truncated Mean", "adjusted truncated Mean", "Median")
  colnames(out) <- c("Estimate", "LB", "UB")
  
  if(!is.finite(tmp.sum$meanWTP)){
    out[1, 2:3] <- -999   # return -999 if the original outcome does not meet the condition for a finite mean WTP estimate 
    output$mWTP <- -999
  }
  
  output$out <- out
  
  class(output) <- "bootCI"
  return(output)
  
}


print.bootCI <- function(x, ...){
  cat("the bootstrap confidence intervals\n")
  printCoefmat(x$out, digits = 5)
  invisible(x$out)
}

# summary.bootCI <- function(object, ...){
#   
#   invisible(object)
# }





