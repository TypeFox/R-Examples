#-->> BEGIN Maximum Likelihood Estimation FOR SUMMARY STATISTICS section

## Generics

setGeneric("cenmle",
           function(obs, censored, groups, ...) standardGeneric("cenmle"))

## Classes

setClass("cenmle", representation("cenreg"))

setClass("cenmle-gaussian", representation("cenmle"))
setClass("cenmle-lognormal", representation("cenmle"))

## Methods

new_cenmle = 
function(x)
{
    if      (is(x, "cenreg-gaussian")) ret = new("cenmle-gaussian", x)
    else if (is(x, "cenreg-lognormal")) ret = new("cenmle-lognormal", x)
    else    stop("Unrecognized distribution") # This shouldn't happen

    return(ret)
}

setMethod("cenmle", 
          signature(obs="Cen", censored="missing", groups="missing"), 
          function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(a~1, list(a=cl[[2]])))
    environment(f) = parent.frame()
    new_cenmle(eval.parent(cenreg(f, ...)))
})

setMethod("cenmle",
          signature(obs="numeric", censored="logical", groups="missing"),
          function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(Cen(a, b)~1, list(a=cl[[2]], b=cl[[3]])))
    environment(f) = parent.frame()
    new_cenmle(eval.parent(cenreg(f, ...)))
})

setMethod("cenmle", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          function(obs, censored, groups, ...)
{
    cl = match.call()
    f = substitute(Cen(a, b)~g, list(a=cl[[2]], b=cl[[3]], g=cl[[4]]))
    f = as.formula(f)
    environment(f) = parent.frame()
    eval.parent(cenreg(f, ...))
})

setMethod("cenmle", 
          signature(obs="numeric", censored="logical", groups="numeric"), 
          function(obs, censored, groups, ...)
{
    cl = match.call()
    f = substitute(Cen(a, b)~g, list(a=cl[[2]], b=cl[[3]], g=cl[[4]]))
    f = as.formula(f)
    environment(f) = parent.frame()
    eval.parent(cenreg(f, ...))
})

setMethod("show", signature(object="cenmle"), function(object)
{
    x = object
    ret = c(x@n, x@n.cen, median(x), mean(x)[1], sd(x)) 
    names(ret) = c("n", "n.cen", "median", "mean", "sd")
    print(ret)
    invisible(ret)
})

setMethod("median", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    as.vector(exp(x@survreg$coef[1]))
})

setMethod("sd", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    coef  = as.vector(x@survreg$coef[1])
    scale = x@survreg$scale[1]

    sqrt(exp(2*coef + scale^2)*(exp(scale^2)-1))
})

setMethod("plot", signature(x="cenmle-lognormal"), 
           function(x, y, xlim=c(-3, 3), ...) 
{
    s = cenfit(x@y, x@ycen)
    plot(x=qnorm(s@survfit$surv), y=log(s@survfit$time), 
         xlab="Normal Quantiles", ylab="log(Value)", xlim=xlim, ...)
    #title(main="Censored Probability Plot")
    abline(x@survreg$coefficients[1], x@survreg$scale)
})

setMethod("plot", signature(x="cenmle-gaussian"), 
           function(x, y, xlim=c(-3, 3), ...) 
{
    s = cenfit(x@y, x@ycen)
    plot(x=qnorm(s@survfit$surv), y=s@survfit$time,
         xlab="Normal Quantiles", ylab="Value", xlim=xlim, ...)
    #title(main="Censored Probability Plot")
    abline(x@survreg$coefficients[1], x@survreg$scale)
})


setMethod("quantile", signature(x="cenmle-lognormal"),
          function(x, probs = NADAprobs, conf.int = FALSE, ...)
{
    q = probs

    int    = as.vector(x@survreg$coef[1])
    scale  = x@survreg$scale

    qhat = exp(int + (qnorm(q) * scale))

    if (!conf.int) 
      {
        ret = qhat
        names(ret) = 
          paste(formatC(100 * q, format = "fg", wid = 1), "%", sep='')
      }
    else
      {
        semean = sqrt(x@survreg$var[1,1])
        cov    = x@survreg$var[2,1] * scale
        varsig = x@survreg$var[2,2] * scale^2

        se = qhat * sqrt(semean^2 + 2*qnorm(q)*cov + (qnorm(q)^2) * varsig)

        # Two-sided conf int
        p = 1-((1-x@conf.int)/2)
        z = qnorm(p)

        w = exp((z*se)/qhat)

        lcl = qhat/w
        ucl = qhat*w

        ret = data.frame(q, qhat, lcl, ucl)
        names(ret) = c("quantile", "value", LCL(x), UCL(x))
      }

    return(ret)
})

setMethod("quantile", signature(x="cenmle-gaussian"),
          function(x, probs = NADAprobs, conf.int = FALSE, ...)
{
    q = probs

    int   = as.vector(x@survreg$coef[1])
    scale = x@survreg$scale

    ncp = qnorm(q)

    qhat = int + (ncp * scale)

    if (!conf.int) ret = qhat
    else
      {
        n  = length(x@survreg$linear.predictors)
        se = sqrt((scale^2)/n)

        # Two-sided conf int
        p = 1-((1-x@conf.int)/2)
        z = qt(p, n-1, abs(ncp))

        lcl = qhat - (z * se)
        ucl = qhat + (z * se)

        ret = data.frame(q, qhat, lcl, ucl)
        names(ret) = c("quantile", "value", LCL(x), UCL(x))
      }

    return(ret)
})

setMethod("mean", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    n     = length(x@survreg$linear.predictors)
    int   = as.vector(x@survreg$coef[1])
    scale = x@survreg$scale

    xbar = as.vector(exp(int + 0.5*(scale)^2))
    se   = sqrt(exp(2*int +scale^2)*(exp(scale^2)-1)/n)

    # Two-sided conf int
    p = 1-((1-x@conf.int)/2)
    gamz = qnorm(p) * sqrt((scale^2/n) + (((0.5)*scale^4)/(n+1)))
    bhat = log(xbar)

    ret = c(xbar, se, exp(bhat - gamz), exp(bhat + gamz))
    names(ret) = c("mean", "se", LCL(x), UCL(x))

    return(ret)
})

setMethod("median", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    as.vector(x@survreg$coef[1])
})

setMethod("sd", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    as.vector(x@survreg$scale)
})

setMethod("mean", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # The mean is the intercept
    int = as.vector(x@survreg$coef[1])

    # The standard error of the intercept/mean
    se = sqrt(x@survreg$var[1,1])

    # Two-sided conf int
    p = 1-((1-x@conf.int)/2)
    gamz = qnorm(p) * se

    ret = c(int, se, (int - gamz), (int + gamz))
    names(ret) = c("mean", "se", LCL(x), UCL(x))

    return(ret)
})

#-->> END Maximum Likelihood Estimation FOR SUMMARY STATISTICS section
