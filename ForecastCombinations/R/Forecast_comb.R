#~~~~~~~~~~~~~Forecast_comb
Forecast_comb <- function(obs, fhat, fhat_new= NULL,
                          Averaging_scheme=c("simple", "ols", "robust", "cls", "variance based", "best")) {
  pckg = c("quantreg", "quadprog")
  temp <- unlist(lapply(pckg, require, character.only=T))
  if (!all(temp==1) ) {
    stop("This function relies on packages \"quadprog\" and \"quantreg\".
         Use ?install.packages if they are not yet installed. \n")
  }
  mat_err <- apply(fhat, 2, function(x) obs - x)
  TT <- NROW(fhat)
  TT_new <- NROW(fhat_new)
  p <- NCOL(fhat)
  pred <- NULL
  weights <- matrix(ncol=p, nrow = TT)
  ## Subroutine needed:
  sq_er <- function(obs, pred) { mean( (obs - pred)^2 )  }

  if(length(Averaging_scheme) != 1) {
    stop("Pick only one of the following:
         c(\"simple\", \"ols\", \"robust\", \"variance based\", \"cls\", \"best\")")
  }

  ## Different forecast averaging schemes
  ## simple
  if(Averaging_scheme== "simple") {
    pred <- apply(fhat, 1, mean)
    weights <- matrix( 1/p, nrow = 1, ncol = p)
    if (!is.null(fhat_new)) { pred_new <- apply(fhat_new, 1, mean) }
    ## OLS weights
  } else if (Averaging_scheme== "ols") {
    weights <- lm(obs ~ fhat)$coef
    pred <- t(weights %*% t(cbind(rep(1, TT), fhat)))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(cbind(rep(1, TT_new), fhat_new))) }
    ## Robust weights
  } else if (Averaging_scheme== "robust") {
    weights <- rq(obs ~ fhat)$coef
    pred <- t(weights %*% t(cbind(rep(1,TT), fhat)))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(cbind(rep(1, TT_new), fhat_new))) }
    ## Based on the variance of the errors. Inverse of MSE
  } else if (Averaging_scheme=="variance based") {
    temp = apply(mat_err^2, 2, mean)/sum(apply(mat_err^2, 2, mean))
    weights <- (1/temp)/sum(1/temp)
    pred <- t(weights %*% t(fhat))
    if ( !is.null(fhat_new) ) { pred_new <- t(weights %*% t(fhat_new)) }
    ## Using constraint least squares
  } else if (Averaging_scheme== "cls") {
    ## Subroutine needed:
    cls1 = function(y, predictions){
      Rinv <- solve(chol(t(predictions) %*% predictions))
      C <- cbind(rep(1, NCOL(predictions)), diag(NCOL(predictions)))
      b = c(1, rep(0, NCOL(predictions)))
      d = t(y) %*% predictions
      qp1 = solve.QP(Dmat= Rinv, factorized= TRUE, dvec= d, Amat= C, bvec = b, meq = 1)
      weights = qp1$sol
      yhatcls = t(weights %*% t(predictions))
      list(yhat= yhatcls, weights= weights)
    }
    ##
    weights <- cls1(obs, fhat)$weights
    pred <- t(weights %*% t(fhat))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(fhat_new)) }
    ## Picking the best performing model up to date according to squared loss function
  } else if (Averaging_scheme== "best") {
    temp <- apply(fhat, 2, sq_er, obs= obs)
    weights <- rep(0, p)
    weights[which.min(temp)] <- 1
    pred <- t(weights %*% t(fhat))
    if (!is.null(fhat_new)) { pred_new <- t(weights %*% t(fhat_new)) }
  }
  if (is.null(fhat_new)) { pred_new <- NULL }
  return( list(fitted= pred, pred = pred_new, weights = weights ) )
  }
