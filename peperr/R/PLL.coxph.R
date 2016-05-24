PLL.coxph <- function(object, newdata, newtime, newstatus, complexity, ...){
 if(is.null(attributes(object)$which.genes)){ # required for Graf&Bauer approach
    logplik(x = newdata, time = newtime, status = newstatus, 
        b = object$coefficients)
    } else {
    z <- matrix(apply(attributes(object)$coef[attributes(object)$which.genes]*t(newdata[,attributes(object)$which.genes]), FUN=sum, MARGIN=2), ncol=1)
    logplik(x = z, time = newtime, status = newstatus,   b = object$coefficients)
    }
}

logplik <- function(x, time, status, b, method = c("breslow", "efron"), return.all=FALSE)
  {
    method <- match.arg(method)
    n <- length(time)
    o <- order(status, decreasing=T)
    oo <- o[order(time[o])]
    time <- time[oo]
    status <- status[oo]
    rept <- rep(0, n)
    for (i in 1:n) rept[i] <- sum(time[i:n]==time[i] & status[i:n]==1)
    complete <- which(status==1)
    nnc <- length(complete)
    dmat <- matrix(0, n, nnc)
    for (i in 1:nnc) {
      dmat[time >= time[complete[i]], i] <- 1
      if (method=="efron") {
        if (rept[complete[i]] > 0) {
          tie <- time==time[complete[i]] & status==1
          di <- max(rept[tie])
          dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]])/di
        }
      }
    }
    eta <- x %*% b
    eeta <- exp(eta)
    k <- ncol(eta)
    loglik <- rep(0, k)
    for (i in 1:k) {
      w <- dmat * eeta[oo, i]
      wsum <- apply(w, 2, sum)
      loglik[i] <- sum(eta[oo, i][status==1]) - sum(log(wsum))
    }
    if (return.all)
      return(list(loglik=loglik, w=scale(w, F, wsum), eta=eta, dmat=dmat, oo=oo))
    else return(loglik)
  }

