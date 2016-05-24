summary.npregress<- function(object, criteria="call", ...) {
  y <- object$call$y
  df <- object$df
  r <- object$residuals
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-df))
  sigma2 <- stderr^2
  if (any(criteria=="call")) {
    criteria <- object$call$criterion
    anscrit <- NULL
  } else {
    crit <-c("aic","aicc","gcv","bic","gmdl")
    if (all(!(criteria%in%crit))) stop(paste("criteria are:",crit,"\n"))
    criteria <- criteria[criteria%in%crit]
    anscrit <- NULL
  }
  if (any(criteria=="gcv"))  anscrit <- c(anscrit,log(sigma2)-2*log(1-df/n))
  if (any(criteria=="aic"))  anscrit <- c(anscrit,log(sigma2)+2*df/n)
  if (any(criteria=="aicc"))  anscrit <- c(anscrit,log(sigma2)+1+(2*(df+1))/(n-df-2))
  if (any(criteria=="bic"))  anscrit <- c(anscrit,log(sigma2) + log(n)*(df)/n)
  if (any(criteria=="gmdl")) {
    Sbul <-   n*sigma2/(n-df)
  anscrit <- c(anscrit,log(Sbul)+df/n*log((sum(y^2)-n*sigma2)/(df*Sbul)))
  }
  if ((criteria!="user")&(!is.null(anscrit))) {
  names(anscrit) <- criteria
} else {
  anscrit <- "No Informative Criterion"
  names(anscrit) <- criteria
}
  ans <- list(residuals=r,Std.Error=stderr,Df=df,Resid.Df=n-df,criteria=anscrit,kernel=object$call$kernel,crit4bw=object$call$criterion,bandwidth=object$bandwidth)
  class(ans) <- "summary.npregress"
  ans
}
