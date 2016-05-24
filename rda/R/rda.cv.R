rda.cv <- function(fit, x, y, prior, alpha, delta,
                   nfold=min(table(y), 10),
                   folds=balanced.folds(y), trace=FALSE){
  this.call <- match.call()

  if(missing(fit)){
    stop("An rda object must be supplied.")
  }

  n <- length(y)
  
  if(missing(folds))
    folds <- split(sample(1:n), rep(1:nfold, length = n))
  else
    nfold <- length(folds)

  if(missing(prior)){
    prior <- fit$prior
  }

  if(missing(alpha)){
    alpha <- fit$alpha
  }

  if(missing(delta)){
    delta <- fit$delta
  }

  err <- cv.err <- ngene <- matrix(0, length(alpha), length(delta))
  dn <- list(round(alpha, 3), round(delta, 3))
  names(dn) <- c("alpha", "delta")
  dimnames(err) <- dimnames(cv.err) <- dimnames(ngene) <- dn

  yhat.new <- array(0, dim=c(length(alpha), length(delta), length(y)))
  dn1 <- list(round(alpha, 3), round(delta, 3), NULL)
  names(dn1) <- c("alpha", "delta", "yhat.new")
  dimnames(yhat.new) <- dn1

  for(k in 1:nfold){
    cat("Fold", k, ":")
    index <- folds[[k]]
    tmp <- rda(x=x[, -index], y=y[-index],
               xnew=x[, index], ynew=y[index],
               alpha=alpha, delta=delta, prior=prior, 
               regularization=fit$reg, trace=trace)
    err <- err+tmp$error
    cv.err <- cv.err+tmp$testerror
    ngene <- ngene+tmp$ngene

    yhat.new[, , index] <- tmp$yhat.new
    cat("\n")
  }

  #dimnames(err) <- dimnames(cv.err) <- list(alpha, delta)
  ngene <- round(ngene/length(folds))
  #dimnames(ngene) <- list(alpha, delta)
  
  obj <- list(alpha=alpha, delta=delta, prior=prior,
              nfold=nfold, folds=folds, yhat.new=yhat.new, 
              err=err, cv.err=cv.err, ngene=ngene,
              call=this.call, reg=fit$reg, n=n)
  class(obj) <- "rdacv"
  obj
}
    
