## Obtain fitted values from ssanova objects
fitted.ssanova <- function(object,...)
{
    mf <- object$mf
    if (!is.null(object$random)) mf$random <- I(object$random$z)
    predict(object,mf)
}

## Obtain residuals from ssanova objects
residuals.ssanova <- function(object,...)
{
    y <- model.response(object$mf,"numeric")
    as.numeric(y-fitted.ssanova(object))
}

## Obtain fitted values in working scale from gssanova objects
fitted.gssanova <- function(object,...)
{
    as.numeric(object$eta)
}

## Obtain residuals from gssanova objects
residuals.gssanova <- function(object,type="working",...)
{
    y <- model.response(object$mf,"numeric")
    wt <- model.weights(object$mf)
    offset <- NULL
    if ((object$family=="nbinomial")&(!is.null(object$nu))) y <- cbind(y,object$nu)
    dat <- switch(object$family,
                  binomial=mkdata.binomial(y,object$eta,wt,offset),
                  nbinomial=mkdata.nbinomial(y,object$eta,wt,offset,object$nu),
                  poisson=mkdata.poisson(y,object$eta,wt,offset),
                  inverse.gaussian=mkdata.inverse.gaussian(y,object$eta,wt,offset),
                  Gamma=mkdata.Gamma(y,object$eta,wt,offset),
                  weibull=mkdata.weibull(y,object$eta,wt,offset,list(object$nu,FALSE)),
                  lognorm=mkdata.lognorm(y,object$eta,wt,offset,list(object$nu,FALSE)),
                  loglogis=mkdata.loglogis(y,object$eta,wt,offset,list(object$nu,FALSE)))
    res <- as.numeric(dat$ywk - object$eta)
    if (!is.na(charmatch(type,"deviance"))) {
        dev.resid <- switch(object$family,
                            binomial=dev.resid.binomial(y,object$eta,wt),
                            nbinomial=dev.resid.nbinomial(y,object$eta,wt),
                            poisson=dev.resid.poisson(y,object$eta,wt),
                            inverse.gaussian=dev.resid.inverse.gaussian(y,object$eta,wt),
                            Gamma=dev.resid.Gamma(y,object$eta,wt),
                            weibull=dev.resid.weibull(y,object$eta,wt,object$nu),
                            lognorm=dev.resid.lognorm(y,object$eta,wt,object$nu),
                            loglogis=dev.resid.loglogis(y,object$eta,wt,object$nu))
        res <- sqrt(dev.resid)*sign(res)
    }
    res
}
