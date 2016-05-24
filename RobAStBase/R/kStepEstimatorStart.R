###############################################################################
## k-step estimator.start
###############################################################################


setMethod("kStepEstimator.start", signature(start = "numeric"),
           function(start, nrvalues, ...){
              if(is.list(start)) start <- unlist(start)
              if(length(start)!=nrvalues)
                  stop("dimension of slot 'estimate' of 'start' != dimension of 'Curve'")
              return(start)
              })

setMethod("kStepEstimator.start", signature(start = "Estimate"),
        function(start, nrvalues, ...){
        if(!is.null(untransformed.estimate(start)))
            return(kStepEstimator.start(untransformed.estimate(start), nrvalues))
        if(!.isUnitMatrix(trafo(start)$mat))
            stop("slot 'untransformed.estimate' of 'start' is null although trafo is non-trivial")
        return(kStepEstimator.start(estimate(start),nrvalues))
})

setMethod("kStepEstimator.start", signature(start = "function"),
        function(start, x, nrvalues, na.rm, L2Fam, startList){
           if(na.rm) x <- na.omit(x)
           start0 <- do.call(start, args=c(list(x,L2Fam),startList))
           return(kStepEstimator.start(start0,nrvalues))
})


