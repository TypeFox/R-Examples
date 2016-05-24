rfit.default <- function (formula, data, subset, yhat0 = NULL, 
    scores = Rfit::wscores, symmetric = FALSE, TAU = 'F0',  ...) {

# Below is taken from quantreg (under GPL) #
  call<-match.call()
  mf<-match.call(expand.dots=FALSE)
  m<-match(c('formula','data','subset'),names(mf),0)
  mf<-mf[c(1,m)]
  mf[[1]]<-as.name("model.frame")
  mf<-eval.parent(mf)
#

  x <- model.matrix(attr(mf, "terms"), data = mf)
  if( abs(max(x) - min(x)) < .Machine$double.eps ^ 0.5 ) stop("x cannot only contain an intercept")
  x1 <- as.matrix(x[,colnames(x)!='(Intercept)'])
  x1 <- as.matrix(cbind(rep(1,nrow(x1)),x1))

  y <- model.response(mf)

  qrx <- qr(x1)
  Q<-as.matrix(qr.Q(qrx))
  q1<-Q[,1]
  xq<-as.matrix(Q[,2:qrx$rank])

  if( is.null(yhat0) ) {
    fit0<-suppressWarnings(rq(y~xq-1))
  } else {
    fit0 <- lsfit(xq, yhat0, intercept = FALSE)
  }
#  ord<-order(fit0$resid)

## 20141211: set initial fit to null model if it has lower dispersion
  betahat0 <- fit0$coef
  if( disp(betahat0, xq, y, scores) > disp(rep(0,length(betahat0)), xq, y, scores) ) {
    betahat0 <- rep(0, length(betahat0) )
  }
  ord <- order(y - xq%*%betahat0)
##

  fit <- jaeckel(as.matrix(xq[ord,]), y[ord], betahat0, scores=scores, ...)
  if( fit$convergence != 0 ) {
    fit2 <- jaeckel(as.matrix(xq[ord,]), y[ord], jitter(fit$par), scores=scores, ...)
    if( fit$convergence != 0 ) {
      warning("rfit: Convergence status not zero in jaeckel")
      if( fit2$value < fit$value ) fit <- fit2
    } else {
      fit <- fit2
    }
    rm(fit2)
  }
  rm(ord)
  betahat <- fit$par

  yhat <- xq %*% betahat
  ehat <- y - yhat
  alphahat <- ifelse(symmetric, signedrank(ehat), median(ehat))
  ehat <- ehat - alphahat
  yhat <- yhat+alphahat

  bhat <- lsfit(x,yhat,intercept=FALSE)$coef

  r.gettau <- switch(TAU,
    F0 = gettauF0,
    R = gettau,
    N = function(...) NA
  )

  tauhat <- r.gettau(ehat, ncol(xq), scores, ...)
  if (symmetric) {
    taushat <- tauhat
  } else {
    taushat <- taustar(ehat, qrx$rank)
  }

  res <- list( coefficients = bhat, residuals = ehat, fitted.values = yhat, 
    scores = scores, x = x, y = y, tauhat = tauhat, qrx1=qrx,
    taushat = taushat, symmetric = symmetric, betahat = bhat,disp=fit$value)
  res$call <- call
  class(res) <- list("rfit")
  res

}
