###
### postMode.R
###

##-----------------------------------------------------------------------------
postMode <- function(y, x, priorCoef) {

tau <- as.double(priorCoef@priorPars['tau'])

if (priorCoef@priorDistr %in% c('pMOM','piMOM','peMOM')) {
  f2min <- function(th) {
    phi <- exp(th[length(th)])
    logl <- sum(dnorm(y,
                      x %*% matrix(th[-length(th)], ncol=1),
                      sd=sqrt(phi),
                      log=TRUE))
    logpr <- priorFun(th)
    return(-logl-logpr)
  }

  ## Set up priorFun
  priorFun <- switch(EXPR=priorCoef@priorDistr,
                     'pMOM' = {
                       r <- as.integer(priorCoef@priorPars['r'])
                       function(th) {
                         dmom(matrix(th[-length(th)], nrow=1),
                              tau=tau,
                              phi=exp(th[length(th)]),
                              r=r,
                              logscale=TRUE,
                              penalty="product")
                       }
                     },
                     'piMOM' = {
                       function(th) {
                         dimom(matrix(th[-length(th)], nrow=1),
                               tau=tau,
                               phi=exp(th[length(th)]),
                               logscale=TRUE,
                               penalty="product")
                       }
                     },
                     'peMOM' = {
                       function(th) {
                         demom(matrix(th[-length(th)], nrow=1),
                               tau=tau,
                               phi=exp(th[length(th)]),
                               logscale=TRUE)
                       }
                     },
                     stop('Prior specified in priorDistr not recognized'))

  ## Optimize
  n <- nrow(x)
  p <- ncol(x)
  lm1 <- lm(y~-1+x)
  thini <- c(coef(lm1), log(sum(residuals(lm1)^2/(n-p))))
  opt <- nlminb(thini, objective=f2min)$par
  ans <- list(coef=opt[-length(opt)],sigma2=exp(opt[length(opt)]))
} else if (priorCoef@priorDistr=='zellner') {
  lm1 <- lm(y~-1+x)
  thopt <- coef(lm1)*tau/(1+tau)
  ans <- list(coef=thopt, sigma2= sum((y - x %*% thopt)^2)/length(y))
} else {
  stop('Prior specified in priorDistr not recognized')
}
return(ans)
}

