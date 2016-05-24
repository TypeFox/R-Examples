###
### Calculate the ATE without covariates
###

ATEnocov <- function(Y, Z, data = parent.frame(), match = NULL){

  ## an internal function that checks match and returns diff
  match.check <- function(Y, Z, match) { 
    n <- length(Y)
    if ((n %% 2) != 0)
      stop("pair randomization requires the even number of observations")
    if (length(unique(table(match))) > 1)
      stop("invalid input for `match'")
    if (unique(table(match)) != 2)
      stop("invalid input for `match'")
    umatch <- sort(unique(match))
    diff <- rep(NA, n/2)
    for (i in 1:length(umatch))
      diff[i] <- Y[(Z == 1) & (match == umatch[i])] -
        Y[(Z == 0) & (match == umatch[i])] 
    return(diff)
  }
  
  ## getting the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  match <- eval(call$match, envir = data)

  ## checking data
  if (sum(sort(unique(Z)) == c(0,1)) != 2)
    stop("`Z' should be binary taking the value of 0 or 1")
  if (length(Y) != length(Z))
    stop("`Y' and `Z' have different numbers of observations")
  if (!is.null(match))
    if (length(match) != length(Y))
      stop("`match' and `Y' have different numbers of observations")
    
  ## ATE for unit randomization
  res$ATE.est <- mean(Y[Z==1])-mean(Y[Z==0])
  res <- list(call = call, Y = Y, Z = Z, match = match)
  if (is.null(match)) { # without matching
    res$ATE.var <- var(Y[Z==1])/sum(Z==1)+var(Y[Z==0])/sum(Z==0)
  } else { # with matching
    res$diff <- diff <- match.check(Y, Z, match)
    res$ATE.var <- var(diff)/length(diff)
  }
  class(res) <- "ATEnocov"
  return(res)
}

