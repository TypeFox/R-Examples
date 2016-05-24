
predict.extraTrees <- function( object, newdata, quantile=NULL, allValues=F, probability=F, newtasks=NULL, ... )
{
    if (!inherits(object, "extraTrees")) {
        stop("Object not of class extraTrees")
    }
    et = object
    
    args <- list(...)
    if (length(args) > 0) stop(sprintf("Illegal argument '%s' to predict.extraTrees.", names(args)[1]))
    
    if (et$multitask && is.null(newtasks)) {
        stop("Prediction requires list of newtasks.")
    }
    
    if (!et$multitask && ! is.null(newtasks)) {
        stop("to predict with newtasks extraTrees must be trained with tasks")
    }
    
    if (probability && ! et$factor) {
        stop("probability=TRUE can be only used for factor model (classification).")
    }
    
    ## making sure no NAs: !!! we now support NAs in Java
    ##if ( any(is.na(newdata)) ) stop("Input matrix newdata contains NAs.")
    newDataNA = any(is.na(newdata))
        
    if (ncol(newdata)!=et$ndim) {
        stop( sprintf("newdata(ncol=%d) does not have the same dimensions as the original x (ncol=%d)", ncol(newdata), et$ndim) )
    }
    if (!is.null(quantile)) {
        if (!et$quantile) {
            stop("to predict quantiles, extraTrees must be trained with quantile=T")
        }
        if (quantile[1] < 0.0 || quantile[1] > 1.0) {
            stop("quantile has to be between 0.0 and 1.0.")
        }
        if (allValues) {
            stop("Can't use allValues=T with quantile.")
        }
        ## quantile regression:
        return( .jcall( et$jobject, "[D", "getQuantiles", toJavaMatrix2D(newdata), quantile[1] ) )
    }
    if (allValues || probability) {
        ## returning allValues prediction:
        if (et$multitask) {
            ## multi-task version:
            m = toRMatrix( .jcall( 
                et$jobject,
                "Lorg/extratrees/data/Matrix;", 
                "getAllValuesMT", 
                toJavaMatrix2D(newdata),
                .jarray(as.integer(newtasks-1))
            ))
        } else {
            m = toRMatrix( .jcall( 
              et$jobject, 
              "Lorg/extratrees/data/Matrix;", 
              "getAllValues", 
              toJavaMatrix2D(newdata) 
            ) )
        }
        ## converting NaN to NA
        m[ is.nan(m) ] = NA
        if (!et$factor) {
            ## regression model:
            return(m)
        }
        if (probability) {
          ## following code assumes (correctly) we have at least 2 classes
          counts = t( apply(m+1, 1, tabulate, nbins=length(et$levels)) )
          counts = counts / rowSums(counts)
          colnames(counts) = et$levels
          return(counts)
        }
        
        ## factor model: convert double matrix into data.frame of factors
        lvls = 0:(length(et$levels)-1)
        ## converting NaN to NA
        m = round(m)
        mlist = lapply( 1:ncol(m), function(j) factor(round(m[,j]), levels=lvls, labels=et$levels) )
        ## changing list to data.frame:
        attributes(mlist) <- list(
            row.names=c(NA_integer_,nrow(m)), 
            class="data.frame",
            names=make.names(names(mlist)),
            unique=TRUE
        )
        return(mlist)
    }
    if ( ! et$factor) {
        if (et$multitask) {
            yhat = .jcall( et$jobject, "[D", "getValuesMT", toJavaMatrix2D(newdata), .jarray(as.integer(newtasks-1)) )
        } else {
            yhat = .jcall( et$jobject, "[D", "getValues", toJavaMatrix2D(newdata) )
        }
        yhat[ is.nan(yhat) ] = NA
        return(yhat)
    }
    if (et$multitask) {
        yhat = .jcall( et$jobject, "[I", "getValuesMT", toJavaMatrix2D(newdata), .jarray(as.integer(newtasks-1)) )
    } else {
        yhat = .jcall( et$jobject, "[I", "getValues", toJavaMatrix2D(newdata) )
    }
    yhat[yhat < 0] = NA
    return( factor(yhat, levels=0:(length(et$levels)-1), labels=et$levels ) )
}
