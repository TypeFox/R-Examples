#-->> BEGIN Regression on Order Statistics (ROS) section

## Generics

setGeneric("ros", 
           function(obs, censored, forwardT = "log", reverseT = "exp") 
             standardGeneric("ros"))

## Classes

setClass("ros", "list")

## Methods

# ros -- Regression on Order Statistics (ROS).
# An implementation of ROS for left-censored data (less-thans) 
# containing one to multiple censoring thresholds. See man page.
setMethod("ros", signature(obs="numeric", censored="logical"),
          function(obs, censored, forwardT = "log", reverseT = "exp")
{
    if (is.null(forwardT) || is.null(reverseT)) {
        forwardT = reverseT = "trueT"
    }
    else if (!exists(forwardT)) {
         stop("Can not find Forward Transformation function: ", forwardT, "\n")
    }
    else if (!exists(reverseT)) {
        stop("Can not find Reverse Transformation function: ", reverseT, "\n")
    }

    if ( (length(obs[censored])/length(obs)) > 0.8 ) {
        warning("Input > 80% censored -- Results are tenuous.\n")
    }

    ix = which(obs > max(obs[!censored]))
    if (length(ix)) {
      warning("Dropped censored values that exceed max of uncensored values.")
      obs = obs[-ix]
      censored = censored[-ix]
    } 

    ix = order(obs)
    obs = obs[ix]
    censored = censored[ix]

    pp = hc.ppoints(obs, censored)

    pp.nq = qnorm(pp[!censored])
    obs.transformed = get(forwardT)(obs[!censored])
    hc = lm(obs.transformed ~ pp.nq)

    oldClass(hc) = c("ros", "lm")
    hc$obs      = obs
    hc$modeled  = obs 
    hc$pp       = pp
    hc$censored = censored 
    hc$reverseT = reverseT
    hc$forwardT = forwardT

    hc$modeled[censored] = predict(hc, qnorm(pp[censored]))

    return(hc)
})

# ROS is technically "Linear Regression on Order Statistics" 
# or LROS.  Version 1.0-1 had funtions and objects as "lros".  
# As of Version 1.0-2, I've renamed things to "ros" to keep consistent
# with the literature on this method which calls it just ROS.
# This is provided for backward compatibility.
lros = ros

# Also, since the rest of the package is evolving to use functions
# names like cen*
cenros = ros

#  trueT is provided so that ros() can be used with no transforms.
#  (It is a quick hack -- rather than recoding ros and family)
trueT =
function(x)
{
    return(x)
}

setMethod("show", signature(object="ros"), function(object)
{
    print(object)
})

setMethod("print", signature(x="ros"), function(x, ...)
{
    n       = length(x$modeled)
    n.cen   = length(x$modeled[x$censored])
    median  = median(x)
    mean    = mean(x)
    sd      = sd(x)

    ret = c(n, n.cen, median, mean, sd)
    names(ret) = c("n", "n.cen", "median", "mean", "sd")

    print(ret, ...)
    invisible(ret)
})

setMethod("summary", "ros", function(object, plot=FALSE, ...)
{
    ret = summary.lm(object, ...)
    if (plot)
      {
        oldClass(object) = "lm"
        plot.lm(object, ...)
      }
    return(ret)
})

# S3 Method -- ros conversion to data.frame discards all linear model info
as.data.frame.ros = 
function (x, row.names = NULL, optional = FALSE, ...)
{
    x = list(obs=x$obs, censored=x$censored, pp=x$pp, modeled=x$modeled)
    as.data.frame(x, row.names, optional)
}

setMethod("mean", signature(x="ros"), function(x, ...) mean(x$modeled, ...))

setMethod("sd", signature(x="ros"), function(x, na.rm=FALSE)
{
    sqrt(var(x=x$modeled))
})

setMethod("median", "ros", function(x, na.rm=FALSE) median(x$modeled))


### min and max methods don't work yet cuz I don't know how to eval
##  the expanded dots correctly yet.
#setMethod("min", "ros", function(..., na.rm=FALSE) 
#{
#    x = match.call(expand.dots=TRUE)[[2]]
#    min(eval.parent(x$modeled))
#})
#setMethod("max", "ros", function(..., na.rm=FALSE) 
#{
#    x = match.call(expand.dots=TRUE)[[2]]
#    max(eval.parent(x$modeled))
#})

## Query and prediction functions for ROS objects

setMethod("quantile", signature(x="ros"), function(x, probs=NADAprobs,...)
{
    quantile(x$modeled, probs, ...)
})

# predict method -- given new normal quantiles of plotting positions,
# returns the corresponding modeled values.
setMethod("predict", "ros", function(object, newdata, ...)
{
    predicted = 
      as.vector(predict.lm(object, newdata=data.frame(pp.nq=newdata), ...))

    return(get(object$reverseT)(predicted))
})

# pexceed method -- given new values or concentrations,
# returns the probability of exceedance.
setMethod("pexceed", "ros",
function(object, newdata, conf.int=FALSE, conf.level=0.95, ...)
{
    cen = object$censored
    forwardT = get(object$forwardT)

    pp.nq = qnorm(object$pp[!cen])
    obs.transformed = forwardT(object$obs[!cen])

    n.lm = lm(pp.nq~obs.transformed)

    nd = data.frame(obs.transformed=forwardT(newdata))
    
    ret = NULL
    if (!conf.int) ret = 1 - pnorm(as.vector(predict.lm(n.lm, nd)))
    else
      {    
        ret = 1 - pnorm(predict.lm(n.lm, nd,
                        interval="confidence", level=conf.level))

        lcl = paste(as.character(conf.level), "LCL", sep='')
        ucl = paste(as.character(conf.level), "UCL", sep='')

        ret[,c(2,3)] = ret[,c(3,2)]
        colnames(ret) = c("prob", lcl, ucl)
        rownames(ret) = as.character(newdata)
      }

    return(ret)
})

## Generic plotting functions

setMethod("lines", signature(x="ros"), function (x, ...)
{
    minmax.nq = qnorm(c(min(x$pp), max(x$pp)))
    lines(x=minmax.nq, y=predict(x, minmax.nq), ...)
})

## Broken for the time being -- use lines
#setMethod("abline", signature(a="ros"), 
#          function(a, b, h, v, reg, coef, untf, col, lty, lwd, ...) 
#{
#    minmax.nq = qnorm(c(min(a$pp), max(a$pp)))
#    lines(x=minmax.nq, y=predict(a, minmax.nq), ...)
#})

# Constructs a prob-plot representation of a ROS model
setMethod("plot", signature(x="ros", y="missing"),
          function(x, 
                   plot.censored = FALSE, 
                   lm.line = TRUE, 
                   grid    = TRUE, 
                   ylab    = "Value", 
                   pch     = 16,
                   ... )
{
    ##
    # To do:
    #   Refactor this routine -- it is long and ugly.
    #   Constrain y ticks to factors of 10?
    #   Draw in grid lines at y ticks
    #
    uncen = x$modeled[!x$censored]
    cen   = x$modeled[x$censored]

    pp.uncen.nq = qnorm(x$pp[!x$censored])
    pp.cen.nq   = qnorm(x$pp[x$censored])

    ymin = min(c(uncen, cen))
    ymax = max(c(uncen, cen))
    xmin = min(c(pp.uncen.nq, pp.cen.nq))
    xmax = max(c(pp.uncen.nq, pp.cen.nq))

    if (x$forwardT == "log") {
        plot(y    = uncen, 
             x    = pp.uncen.nq, 
             ylim = c(ymin, ymax), 
             xlim = c(xmin, xmax), 
             ylab = ylab,
             xlab = "Normal Quantiles",
             log  = "y",
             yaxt = "n",
             pch  = 16,
             ...
        )
     }
     else {
        plot(y    = uncen, 
             x    = pp.uncen.nq, 
             ylim = c(ymin, ymax), 
             xlim = c(xmin, xmax), 
             ylab = ylab,
             xlab = "Normal Quantiles",
             yaxt = "n",
             pch  = 16,
             ...
        )
     }

    axis(2, axTicks(2))

    # Draw a line through the predicted xmin and xmax points 
    if (lm.line) lines(x)

    # Draw top "Chance of Exceedance" axis
    labels = c("5","10","25","50","75","90","95")
    atv = qnorm(c(.05, .1, .25, .50, .75, .90, .95))
    axis(3, at=atv, labels=rev(labels), las=2)
    mtext("Percent Chance of Exceedance", side=3, line=3)

    # Plot the synthetic censored points -- if requested
    if (plot.censored) points(y=cen, x=pp.cen.nq)

    # Draw in grid at major divisions  -- if requested
    if (grid)
      {
        abline(v=atv, lty="dotted")
        abline(h=axTicks(2), lty="dotted")
      }
})

#setMethod("boxplot", signature(x="ros"), function(x, log="y", range=0,...)
setMethod("boxplot", signature(x="ros"), function(x, ...)
{
    # use boxplot.default instead
    bx = boxplot(x$modeled, plot=FALSE, ...)

    # The vector stats is the famous five -- with outlier limts.
    # We relpace the quantiles with the ros-modeled ones
    bx$stats[2:4] = as.vector(quantile(x, c(0.25, 0.5, 0.75)))
    bx$stats = as.matrix(bx$stats)

    bxp(bx, ...)
    invisible(bx)
})

## Routines for calculating Helsel-Cohn style plotting positions

# hc.ppoints calculates computes Wiebull-type plotting postions of data
# containing mixed uncensored and left-censored data. See man page.
# If there are no censored values, the plotting postitions are calculated
# using the standard ppoints() function.
hc.ppoints = 
function(obs, censored)
{    
    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    pp = numeric(length(obs))

    if (!any(censored)) pp = ppoints(obs)
    else
      {
        #cn = cohn(obs, censored)
        #pp[!censored] = hc.ppoints.uncen(cn=cn)
        #pp[censored]  = hc.ppoints.cen(cn=cn)
        pp[!censored] = hc.ppoints.uncen(obs, censored)
        pp[censored]  = hc.ppoints.cen(obs, censored)
      }

    return(pp)
}

# cohn Calculates "Cohn" Numbers -- quantities described by
# Helsel and Cohn's (1988) reformulation of a prob-plotting formula
# described by Hirsch and Stedinger (1987).
#
# The Cohn Numbers are:
# A_j   = the number of uncensored obs above the jth threshold.
# B_j   = the number of observations (cen & uncen) below the jth threshold.
# C_j   = the number of censored observations at the jth threshold
# P_j   = the probability of exceeding the jth threshold

cohn =
function(obs, censored)
{
     uncen = obs[!censored]
     cen   = obs[censored]

     A = B = C = P = numeric()

     limit = sort(unique(cen))

     a = length(uncen[uncen < limit[1]])
     if (a > 0) { limit = c(0, limit) }

     i = length(limit)

     A[i] = length(uncen[ uncen >= limit[i] ])
     B[i] = length(obs[obs <= limit[i]])-length(uncen[uncen==limit[i]])
     C[i] = length(cen[ cen == limit[i] ])
     P[i] = A[i]/(A[i] + B[i])

     i = i - 1
     while (i > 0)
       {
         A[i] = length(uncen[ uncen >= limit[i] & uncen < limit[i + 1] ])
         B[i] = length(obs[obs <= limit[i]])-length(uncen[uncen==limit[i]])
         C[i] = length(cen[cen == limit[i]])
         P[i] = P[i + 1] + ((A[i]/(A[i] + B[i])) * (1 - P[i + 1]))

         i = i - 1
       }

     return(list(A=A, B=B, C=C, P=P, limit=limit))
}

# hc.ppoints.uncen calculates plotting postions for uncensored data.
hc.ppoints.uncen =
function(obs, censored, cn)
{
    #if (missing(cn)) { cn = cohn(obs, censored) }
    cn = cohn(obs, censored)

    nonzero = (cn$A != 0)
    A     = cn$A[nonzero]
    B     = cn$B[nonzero]
    P     = cn$P[nonzero]
    limit = cn$limit[nonzero]

    pp = numeric()
    n = length(limit)
    for (i in 1:n)
      {
        R = 1:A[i] 

        k = P[i+1]
        if (is.na(k)) k = 0

        for (r in 1:length(R))
          {
            pp = c(pp, (1 - P[i]) + ((P[i] - k) * R[r])/(A[i] + 1))
          }
      }
    return(pp)
}

# hc.ppoints.cen calculates plotting postions for censored data.
hc.ppoints.cen =
function(obs, censored, cn)
{    
    if (missing(cn)) { cn = cohn(obs, censored) }

    C     = cn$C
    P     = cn$P
    limit = cn$limit

    if (P[1] == 1) 
      {
        C     = C[-1]
        P     = P[-1]
        limit = limit[-1]
      }
    
    pp = numeric()
    for (i in 1:length(limit)) 
      {
        c.i = C[i]
        for (r in 1:c.i) 
          {
            pp = c(pp, (1 - P[i]) * r/(c.i + 1))
          }
      }
    return(pp)
}

#-->> END Regression on Order Statistics (ROS) section
