## METHODS FOR EXTRACTING INFORMATION FROM THE nplr CLASS
setMethod("getX", "nplr", function(object) return(object@x))
setMethod("getY", "nplr", function(object) return(object@y))
setMethod("getW", "nplr", function(object) return(object@w))
setMethod("getFitValues", "nplr", function(object) return(object@yFit))
setMethod("getXcurve", "nplr", function(object) return(object@xCurve))
setMethod("getYcurve", "nplr", function(object) return(object@yCurve))
setMethod("getInflexion", "nplr", function(object) return(object@inflPoint))
setMethod("getPar", "nplr", function(object){
    return(list(npar=object@npars, params=object@pars))
    })
setMethod('getGoodness', 'nplr', function(object) return(object@goodness))
setMethod('getStdErr', 'nplr', function(object) return(object@stdErr))
setMethod("getNlmErr", "nplr", function(object) return(object@nlmErr))
setMethod("getAUC", "nplr", function(object) return(object@AUC))


## MAIN nplr FUNCION
nplr <- function(x, y, useLog=TRUE, LPweight=0.25,
                npars="all", method=c("res", "sdw", "gw"),
                silent=FALSE){

    if(length(x)!=length(y))
        stop("x and y lengths differ.")
  
    if(is.numeric(npars) & (npars<2 | npars>5))
        stop("\n'npars' must be in [2, 5], or 'all'!")
  
    method <- match.arg(method)

    repTable <- table(x)
    maxrep <- max(repTable, na.rm=TRUE)
    minrep <- min(repTable, na.rm=TRUE)
    if(method=="sdw"){
        if(maxrep<2){
            method <- "res"
            if(!silent){
                warning("\nNone of the x-values seem to be replicated.
                    The 'sdw' method has been replaced by 'res'.",
                    call.=FALSE, immediate.=TRUE)
                message()
            }

        } else if(minrep<2){
            if(!silent){
                warning("\nOne (or more) points have no replicates.
                    The 'sdw' method may not be appropriate.",
                    call.=FALSE, immediate.=TRUE)
                message()
            }
            }
    }
  
    if(method=="gw" & any(y<0))
        if(!silent){
            warning("\nBecause of one (or more) y negative values,
                the 'gw' method may not be appropriate.",
                call.=FALSE, immediate.=TRUE)
            message()
        }
    if(any(is.na(x) | is.na(y))){
        NAs <- union(which(is.na(x)), which(is.na(y)))
        x <- x[-NAs]
        y <- y[-NAs]
        if(!silent){
            warning(call.=FALSE,
                sprintf("%s point(s) has(ve) been removed for missingness.",
                    length(NAs)),
                immediate.=TRUE)
            message()
                }
    }
    y <- y[order(x)]
    x <- sort(x)
  
    pp <- sum(y<0 | y>1)/length(y)
    if(pp > .2 & !silent){
        warningtext <- "% of your y values fall outside the range [0, 1]"
        warning(call.=FALSE,
            sprintf("%s%s", round(pp*100, 2), warningtext),
            immediate.=TRUE)
        message("\t- any results output may not be representative.")
        message("\t- be sure you are using y-values as proportions.")
        message()
    }
  
    if(useLog) x <- log10(x)
    object <- new("nplr", x=x, y=y, useLog=useLog, LPweight=LPweight)
    object@call <- match.call()
  
    .weight <- .chooseWeight(method)

    if(npars=="all"){
        testAll <- .testAll(.sce, x, y, .weight, LPweight, silent)
        npars <- testAll$npars
        if(!silent){
            msg <- sprintf("The %s-parameters model showed better performance",
                format(npars))
            message(msg)
        }
    }
  
    nPL <- .chooseModel(npars)
    inits <- .initPars(x, y, npars)
    options(warn = -1)
    best <- nlm(f=.sce, p=inits, x=x, yobs=y, .weight, LPweight, nPL)
    options(warn = 0)
    if(best$iterations==0)
        stop("'nlm' failed to estimate parameters.\n")
  
    # Best estimates
    bottom <- best$estimate[1]
    top <- best$estimate[2]
    xmid<-best$estimate[3]
    scal <- best$estimate[4]
    s <- best$estimate[5]
  
    # Estimating values
    newX <- seq(min(x), max(x), length=200)
    newY <- nPL(bottom, top, xmid, scal, s, newX)
    yFit <- nPL(bottom, top, xmid, scal, s, x)
    if(length(unique(signif(yFit, 5)))==1)
        stop("nplr failed and returned constant fitted values.
            Your data may not be appropriate for such model.")

    perf <- .getPerf(y, yFit)
    
    # Inflexion point coordinates
    pars <- cbind.data.frame(bottom=bottom, top=top, xmid=xmid, scal=scal, s=s)
    infl <- .inflPoint(pars)
  
    object@w <- (y - yFit)^2
#    object@nlmErr <- ifelse(npars=="all", testAll$err, NA)
    object@npars <- npars
    object@pars <- pars
    object@yFit <- yFit
    object@xCurve <- newX
    object@yCurve <- newY
    object@inflPoint <- infl
    object@goodness <- perf$goodness
    object@stdErr <- perf$stdErr
    object@AUC <- data.frame(trapezoid = .AUC(newX, newY), Simpson = .Simpson(newX, newY))
  
    return(object)
}
