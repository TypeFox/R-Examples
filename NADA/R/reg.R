#-->> BEGIN Maximum Likelihood Estimation REGRESSION section

## Generics

setGeneric("cenreg",
  function(obs, censored, groups, ...) standardGeneric("cenreg"))

## Classes

setOldClass("survreg")

setClass("cenreg", 
         representation(
          y="numeric", ycen="logical", 
          n="numeric", n.cen="numeric", conf.int="numeric", survreg="survreg"))

setClass("cenreg-gaussian", representation("cenreg"))
setClass("cenreg-lognormal", representation("cenreg"))

setClass("summary.cenreg", representation("list"))

## Methods

setMethod("LCL", 
          signature(x="cenreg"), function(x) paste(x@conf.int, "LCL", sep=''))

setMethod("UCL", 
          signature(x="cenreg"), function(x) paste(x@conf.int, "UCL", sep=''))

# The generic formula method that is called from below
cenreg.default =
function(obs, censored, groups, dist, conf.int=0.95, ...)
{
    dist = ifelse(missing(dist), "lognormal", dist)
    ret = switch(dist,
                 gaussian  = new_cenreg_gaussian(obs, ...),
                 lognormal = new_cenreg_lognormal(obs, ...),
                 survreg(asSurv(obs), dist=dist, ...)
    )

    ret@y    = eval(obs[[2]][[2]], environment(obs))
    ret@ycen = eval(obs[[2]][[3]], environment(obs))

    ret@n = length(ret@y)
    ret@n.cen = length(ret@y[ret@ycen])

    ret@conf.int = conf.int

    #browser()
    ret@survreg$call = sys.call(sys.parent())

    return(ret)
}

setMethod("cenreg",
          signature(obs="formula", censored="missing", groups="missing"),
          cenreg.default)

setMethod("cenreg", 
          signature(obs="numeric", censored="logical", groups="numeric"), 
          cencen.vectors.groups)

setMethod("cenreg", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          cencen.vectors.groups)

setMethod("summary", signature(object="cenreg"), function(object, ...)
{
    s = unclass(summary(object@survreg, ...))
    s$r = rloglik(s)
    return (new("summary.cenreg", s))
})

setMethod("show", signature(object="cenreg"), function(object)
{
    ret = summary(object)
    print(ret)
    invisible(ret)
})

setMethod("plot", signature(x="cenreg"), 
           function(x, y, xlab="Normal Quantiles", ylab="Residuals", 
                    xlim=c(-xmax, xmax), ylim=c(-ymax, ymax), ...) 
{
    res = as.vector(residuals(x, type="deviance"))
    res.km = cenfit(res, x@ycen)

    surv = res.km@survfit$surv
    time = res.km@survfit$time

    nq = qnorm(surv)
    xmax = max(nq)
    ymax = max(time)
    plot(nq, time, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, ...)
})

setMethod("predict", signature(object="cenreg-lognormal"), 
          function(object, newdata, q=NADAprobs, ...)
{

    p      = c(1, newdata)
    coefs  = as.vector(object@survreg$coefficients)

    if (length(p) != length(coefs))
      {
        stop("input predictors != number of model parameters")
      }

    df     = object@survreg$df
    var    = object@survreg$var[-df, -df]
    scale  = object@survreg$scale
    varsig = object@survreg$var[df, df] * scale^2
    cov    = p %*% (object@survreg$var[df,][-df] * scale)

    varmean = t(p) %*% var %*% p

    qhat = exp(sum(p * coefs) + (qnorm(q) * scale))

    se = qhat * sqrt(varmean + (2*qnorm(q)*cov) + (qnorm(q)^2) * varsig)

    ci =  1-((1-object@conf.int)/2)
    z = qnorm(ci)

    w = exp((z*se)/qhat)

    lcl = qhat/w
    ucl = qhat*w

    ret = data.frame(q, qhat, se, lcl, ucl)
    colnames(ret) = c("quantile", "value", "se", LCL(object), UCL(object))

    return(ret)
})

setMethod("residuals", signature(object="cenreg"), function(object, ...)
{
    residuals(object@survreg, ...)
})

setMethod("coef", signature(object="cenreg"), function(object, ...)
{
    coef(object@survreg, ...)
})

setMethod("cor", signature(x="cenreg"), function(x, y, use, method)
{
    # Calls private supporting function (below)
    rloglik(summary(x))
})

# This is summary.survreg from the survival package --
# for now we hack it to do what we want. 
setMethod("show", signature(object="summary.cenreg"), function(object)
{
    x = object
    digits = max(options()$digits - 4, 3)
    n <- x$n

    print(x$table, digits = digits)
    if (nrow(x$var)==length(x$coefficients)) 
      {
	    cat("\nScale fixed at", format(x$scale, digits=digits),"\n") 
      }
    else 
      {
	    cat ("\nScale =", format(x$scale, digits=digits), "\n")
      }

    cat("\n", x$parms, "\n", sep='')
    df  <- sum(x$df) - x$idf   # The sum is for penalized models
    cat("Loglik(model)=", format(round(x$loglik[2],1)),
	"  Loglik(intercept only)=", format(round(x$loglik[1],1)), "\n")

    cat("Loglik-r: ", x$r, "\n");

    if (df <= 0) cat ("\n")
    else
      {
	    cat("\nChisq=", format(round(x$chi,2)), "on", round(df,1),
		"degrees of freedom, p=", 
		format(signif(1-pchisq(x$chi, df),2)), "\n")
      }

    if (x$robust) cat("(Loglikelihood assumes independent observations)\n")

    cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)), "\n")

    if (!length(x$na.action)) cat("n =", x$n, "\n")
	else cat("n =", x$n, " (", naprint(x$na.action), ")\n", sep="")

    if(!is.null(x$correl)) 
      {
        p <- dim(x$correl)[2]
        if(p > 1) 
          {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(x$correl)
            x$correl[ll] <- format(round(x$correl[ll], digits=digits))
            x$correl[!ll] <- ""
            print(x$correl[-1,  - p, drop = FALSE], quote = FALSE)
          }
      }
    cat("\n")
    invisible(NULL)
})

## Supporting Functions -- private

# Compute the log-likelihood correlation coef (r) from a cenreg summary obj.
# Eq 11.4 pg 187 of Dennis' book
rloglik =
function(x)
{
    n = x$n
    G = -2 * diff(x$loglik)
    sqrt(1 - exp(G/n))
}

# cenreg for lognormal distributions

new_cenreg_lognormal =
function(formula, dist, ...)
{
    s = survreg=survreg(asSurv(formula), dist="lognormal", ...)
    new("cenreg-lognormal", survreg=s)
}

# cenreg for gaussian, or normal, distributions

# If a normal distribution is assumed the input data must be expressed
# as an interval between zero and the DL.  They cannot simply be stated as
# 'left' censored, because that allows some probability of going below 0.
# Since environmental/analytical data are _usually_ not negative,
# estimates will be biased low and wrong.  So with the normal option
# and left censoring, internally we must use interval censoring.  The end of
# the interval are the detected values.  The start of the interval will
# have identical numbers in it for the detects, and a 0 for the 
# nondetects (a simple trick is: start = obs - obs * censored).

new_cenreg_gaussian =
function(formula, ...)
{
    end      = eval(formula[[2]][[2]], environment(formula))
    censored = eval(formula[[2]][[3]], environment(formula))

    start = end - (end * censored)
    groups = formula[[3]]
    environment(groups) = environment(formula)

    f = substitute(Surv(start, end, type="interval2")~g, list(g=groups))
    f = as.formula(f)

    cl = call("survreg", f, dist="gaussian")
    s = eval(substitute(cl), list(start=start, end=end), environment(formula))
      
    new("cenreg-gaussian", survreg=s)
}

#-->> END Maximum Likelihood Estimation REGRESSION section
