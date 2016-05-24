performance <- function(prediction.obj, measure,
                        x.measure="cutoff", ...) {

    ## define the needed environments
    envir.list <- .define.environments()
    long.unit.names <- envir.list$long.unit.names
    function.names <- envir.list$function.names
    obligatory.x.axis <- envir.list$obligatory.x.axis
    optional.arguments <- envir.list$optional.arguments
    default.values <- envir.list$default.values
    
    ## abort in case of misuse
    if (class(prediction.obj) != 'prediction' ||
        !exists(measure, where=long.unit.names, inherits=FALSE) ||
        !exists(x.measure, where=long.unit.names, inherits=FALSE)) {
      stop(paste("Wrong argument types: First argument must be of type",
                 "'prediction'; second and optional third argument must",
                 "be available performance measures!"))
    }
    
    ## abort, if attempt is made to use a measure that has an obligatory
    ## x.axis as the x.measure (cannot be combined)
    if (exists( x.measure, where=obligatory.x.axis, inherits=FALSE )) {
        message <- paste("The performance measure",
                         x.measure,
                         "can only be used as 'measure', because it has",
                         "the following obligatory 'x.measure':\n",
                         get( x.measure, envir=obligatory.x.axis))
        stop(message)
    }

    ## if measure is a performance measure with obligatory x.axis, then
    ## enforce this axis:
    if (exists( measure, where=obligatory.x.axis, inherits=FALSE )) {
        x.measure <- get( measure, envir=obligatory.x.axis )
    }

    if (x.measure == "cutoff" ||
        exists( measure, where=obligatory.x.axis, inherits=FALSE )) {

        ## fetch from '...' any optional arguments for the performance
        ## measure at hand that are given, otherwise fill up the default values
        optional.args <- list(...)
        argnames <- c()
        if ( exists( measure, where=optional.arguments, inherits=FALSE )) {
            argnames <- get( measure, envir=optional.arguments )
            default.arglist <- list()
            for (i in 1:length(argnames)) {
                default.arglist <- c(default.arglist,
                                     get(paste(measure,":",argnames[i],sep=""),
                                         envir=default.values, inherits=FALSE))
            }
            names(default.arglist) <- argnames

            for (i in 1:length(argnames)) {
                templist <- list(optional.args,
                                 default.arglist[[i]])
                names(templist) <- c('arglist', argnames[i])
                
                optional.args <- do.call('.farg', templist)
            }
        }
        optional.args <- .select.args( optional.args, argnames )
        
        ## determine function name
        function.name <- get( measure, envir=function.names )

        ## for each x-validation run, compute the requested performance measure
        x.values <- list()
        y.values <- list()
        for (i in 1:length( prediction.obj@predictions )) {
            argumentlist <- .sarg(optional.args,
                                  predictions= prediction.obj@predictions[[i]],
                                  labels= prediction.obj@labels[[i]],
                                  cutoffs= prediction.obj@cutoffs[[i]],
                                  fp= prediction.obj@fp[[i]],
                                  tp= prediction.obj@tp[[i]],
                                  fn= prediction.obj@fn[[i]],
                                  tn= prediction.obj@tn[[i]],
                                  n.pos= prediction.obj@n.pos[[i]],
                                  n.neg= prediction.obj@n.neg[[i]],
                                  n.pos.pred= prediction.obj@n.pos.pred[[i]],
                                  n.neg.pred= prediction.obj@n.neg.pred[[i]])

            ans <- do.call( function.name, argumentlist )

            if (!is.null(ans[[1]])) x.values <- c( x.values, list( ans[[1]] ))
            y.values <- c( y.values, list( ans[[2]] ))
        }

        if (! (length(x.values)==0 || length(x.values)==length(y.values)) ) {
            stop("Consistency error.")
        }
        
        ## create a new performance object
        return( new("performance",
                    x.name       = get( x.measure, envir=long.unit.names ),
                    y.name       = get( measure, envir=long.unit.names ),
                    alpha.name   = "none",
                    x.values     = x.values,
                    y.values     = y.values,
                    alpha.values = list() ))
    } else {
        perf.obj.1 <- performance( prediction.obj, measure=x.measure, ... )
        perf.obj.2 <- performance( prediction.obj, measure=measure, ... )
        return( .combine.performance.objects( perf.obj.1, perf.obj.2 ) )
    }
}

.combine.performance.objects <- function( p.obj.1, p.obj.2 ) {
    ## some checks for misusage (in any way, this function is
    ## only for internal use)
    if ( p.obj.1@x.name != p.obj.2@x.name ) {
        stop("Error: Objects need to have identical x axis.")
    }
    if ( p.obj.1@alpha.name != "none" || p.obj.2@alpha.name != "none") {
        stop("Error: At least one of the two objects has already been merged.")
    }
    if (length(p.obj.1@x.values) != length(p.obj.2@x.values)) {
        stop(paste("Only performance objects with identical number of",
                   "cross-validation runs can be combined."))
    }

    x.values <- list()
    x.name <- p.obj.1@y.name
    y.values <- list()
    y.name <- p.obj.2@y.name
    alpha.values <- list()
    alpha.name <- p.obj.1@x.name

    for (i in 1:length( p.obj.1@x.values )) {
        x.values.1 <- p.obj.1@x.values[[i]]
        y.values.1 <- p.obj.1@y.values[[i]]
        x.values.2 <- p.obj.2@x.values[[i]]
        y.values.2 <- p.obj.2@y.values[[i]]

        ## cutoffs of combined object = merged cutoffs of simple objects
        cutoffs <- sort( unique( c(x.values.1, x.values.2)), decreasing=TRUE )

        ## calculate y.values at cutoffs using step function
        y.values.int.1 <- approxfun(x.values.1, y.values.1,
                                    method="constant",f=1,rule=2)(cutoffs)
        y.values.int.2 <- approxfun(x.values.2, y.values.2,
                                    method="constant",f=1,rule=2)(cutoffs)

        ## 'approxfun' ignores NA and NaN
        objs <- list( y.values.int.1, y.values.int.2)
        objs.x <- list( x.values.1, x.values.2 )
        na.cutoffs.1.bool <- is.na( y.values.1) & !is.nan( y.values.1 )
        nan.cutoffs.1.bool <- is.nan( y.values.1)
        na.cutoffs.2.bool <- is.na( y.values.2) & !is.nan( y.values.2 )
        nan.cutoffs.2.bool <- is.nan( y.values.2)
        bools <- list(na.cutoffs.1.bool, nan.cutoffs.1.bool,
                      na.cutoffs.2.bool, nan.cutoffs.2.bool)
        values <- c(NA,NaN,NA,NaN)
        
        for (j in 1:4) {
            for (k in which(bools[[j]])) {
                interval.max <- objs.x[[ ceiling(j/2) ]][k]
                interval.min <- -Inf
                if (k < length(objs.x[[ ceiling(j/2) ]])) {
                    interval.min <- objs.x[[ ceiling(j/2) ]][k+1]
                }
                objs[[ ceiling(j/2) ]][cutoffs <= interval.max &
                                       cutoffs > interval.min ] <- values[j]
            }
        }

        alpha.values <- c(alpha.values, list(cutoffs))
        x.values <- c(x.values, list(objs[[1]]))
        y.values <- c(y.values, list(objs[[2]]))
    }
    
    return( new("performance",
                x.name=x.name, y.name=y.name,
                alpha.name=alpha.name, x.values=x.values,
                y.values=y.values, alpha.values=alpha.values))
}

.define.environments <- function() {
    ## There are five environments: long.unit.names, function.names,
    ## obligatory.x.axis, optional.arguments, default.values
    
    ## Define long names corresponding to the measure abbreviations.
    long.unit.names <- new.env()
    assign("none","None", envir=long.unit.names)
    assign("cutoff", "Cutoff", envir=long.unit.names)
    assign("acc", "Accuracy", envir=long.unit.names)
    assign("err", "Error Rate", envir=long.unit.names)
    assign("fpr", "False positive rate", envir=long.unit.names)
    assign("tpr", "True positive rate", envir=long.unit.names)
    assign("rec", "Recall", envir=long.unit.names)
    assign("sens", "Sensitivity", envir=long.unit.names)
    assign("fnr", "False negative rate", envir=long.unit.names)
    assign("tnr", "True negative rate", envir=long.unit.names)
    assign("spec", "Specificity", envir=long.unit.names)
    assign("ppv", "Positive predictive value", envir=long.unit.names)
    assign("prec", "Precision", envir=long.unit.names)
    assign("npv", "Negative predictive value", envir=long.unit.names)
    assign("fall", "Fallout", envir=long.unit.names)
    assign("miss", "Miss", envir=long.unit.names)
    assign("pcfall", "Prediction-conditioned fallout", envir=long.unit.names)
    assign("pcmiss", "Prediction-conditioned miss", envir=long.unit.names)
    assign("rpp", "Rate of positive predictions", envir=long.unit.names)
    assign("rnp", "Rate of negative predictions", envir=long.unit.names)
    assign("auc","Area under the ROC curve", envir=long.unit.names)
    assign("cal", "Calibration error", envir=long.unit.names)
    assign("mwp", "Median window position", envir=long.unit.names)
    assign("prbe","Precision/recall break-even point", envir=long.unit.names)
    assign("rch", "ROC convex hull", envir=long.unit.names)
    assign("mxe", "Mean cross-entropy", envir=long.unit.names)
    assign("rmse","Root-mean-square error", envir=long.unit.names)
    assign("phi", "Phi correlation coefficient", envir=long.unit.names)
    assign("mat","Matthews correlation coefficient", envir=long.unit.names)
    assign("mi", "Mutual information", envir=long.unit.names)
    assign("chisq", "Chi-square test statistic", envir=long.unit.names)
    assign("odds","Odds ratio", envir=long.unit.names)
    assign("lift", "Lift value", envir=long.unit.names)
    assign("f","Precision-Recall F measure", envir=long.unit.names)
    assign("sar", "SAR", envir=long.unit.names)
    assign("ecost", "Expected cost", envir=long.unit.names)
    assign("cost", "Explicit cost", envir=long.unit.names)

    ## Define function names corresponding to the measure abbreviations.
    function.names <- new.env()
    assign("acc", ".performance.accuracy", envir=function.names)
    assign("err", ".performance.error.rate", envir=function.names)
    assign("fpr", ".performance.false.positive.rate", envir=function.names)
    assign("tpr", ".performance.true.positive.rate", envir=function.names)
    assign("rec", ".performance.true.positive.rate", envir=function.names)
    assign("sens", ".performance.true.positive.rate", envir=function.names)
    assign("fnr", ".performance.false.negative.rate", envir=function.names)
    assign("tnr", ".performance.true.negative.rate", envir=function.names)
    assign("spec", ".performance.true.negative.rate", envir=function.names)
    assign("ppv", ".performance.positive.predictive.value",
           envir=function.names)
    assign("prec", ".performance.positive.predictive.value",
           envir=function.names)
    assign("npv", ".performance.negative.predictive.value",
           envir=function.names)
    assign("fall", ".performance.false.positive.rate", envir=function.names)
    assign("miss", ".performance.false.negative.rate", envir=function.names)
    assign("pcfall", ".performance.prediction.conditioned.fallout",
           envir=function.names)
    assign("pcmiss", ".performance.prediction.conditioned.miss",
           envir=function.names)
    assign("rpp", ".performance.rate.of.positive.predictions",
           envir=function.names)
    assign("rnp", ".performance.rate.of.negative.predictions",
           envir=function.names)
    assign("auc", ".performance.auc", envir=function.names)
    assign("cal", ".performance.calibration.error", envir=function.names)
    assign("prbe", ".performance.precision.recall.break.even.point",
           envir=function.names)
    assign("rch", ".performance.rocconvexhull", envir=function.names)
    assign("mxe", ".performance.mean.cross.entropy", envir=function.names)
    assign("rmse", ".performance.root.mean.squared.error",
           envir=function.names)
    assign("phi", ".performance.phi", envir=function.names)
    assign("mat", ".performance.phi", envir=function.names)
    assign("mi", ".performance.mutual.information", envir=function.names)
    assign("chisq", ".performance.chisq", envir=function.names)
    assign("odds", ".performance.odds.ratio", envir=function.names)
    assign("lift", ".performance.lift", envir=function.names)
    assign("f", ".performance.f", envir=function.names)
    assign("sar", ".performance.sar", envir=function.names)
    assign("ecost", ".performance.expected.cost", envir=function.names)
    assign("cost", ".performance.cost", envir=function.names)

    ## If a measure comes along with an obligatory x axis (including "none"),
    ## list it here.
    obligatory.x.axis <- new.env()
    assign("mxe", "none", envir=obligatory.x.axis)
    assign("rmse", "none", envir=obligatory.x.axis)
    assign("prbe", "none", envir=obligatory.x.axis)
    assign("auc", "none", envir=obligatory.x.axis)
    assign("rch","none", envir=obligatory.x.axis)
    ## ecost requires probability cost function as x axis, which is handled
    ## implicitly, not as an explicit performance measure.
    assign("ecost","none", envir=obligatory.x.axis)  
    
    ## If a measure has optional arguments, list the names of the
    ## arguments here.
    optional.arguments <- new.env()
    assign("cal", "window.size", envir=optional.arguments)
    assign("f", "alpha", envir=optional.arguments)
    assign("cost", c("cost.fp", "cost.fn"), envir=optional.arguments)
    assign("auc", "fpr.stop", envir=optional.arguments)
        
    ## If a measure has additional arguments, list the default values
    ## for them here. Naming convention: e.g. "cal" has an optional
    ## argument "window.size" the key to use here is "cal:window.size"
    ## (colon as separator)
    default.values <- new.env()
    assign("cal:window.size", 100, envir=default.values)
    assign("f:alpha", 0.5, envir=default.values)
    assign("cost:cost.fp", 1, envir=default.values)
    assign("cost:cost.fn", 1, envir=default.values)
    assign("auc:fpr.stop", 1, envir=default.values) 
    
    list(long.unit.names=long.unit.names, function.names=function.names,
         obligatory.x.axis=obligatory.x.axis,
         optional.arguments=optional.arguments,
         default.values=default.values)
}
