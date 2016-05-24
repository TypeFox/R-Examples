

setMethod(f = "classif",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(Classif(object)) } )

setMethod(f = "coef",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(Coef(object)) } )

setMethod(f = "estimator",
          signature = c(x="DynTxRegime"),
          definition = function(x, ...){ return(Est(x,...)) } )

setMethod(f = "fitObject",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(FitObject(object, ...)) } )

setMethod(f = "FM",
          signature = c(object="IQLearnSS"),
          definition = function(object, ...){ return(FittedMain(object, ...)) } )

setMethod(f = "FC",
          signature = c(object="IQLearnSS"),
          definition = function(object, ...){ return(FittedCont(object, ...)) } )

setMethod(f = "FM",
          signature = c(object="IQLearnFS"),
          definition = function(object, ...){ return(FittedMain(object, ...)) } )

setMethod(f = "FC",
          signature = c(object="IQLearnFS"),
          definition = function(object, ...){ return(FittedCont(object, ...)) } )

setMethod(f = "FM",
          signature = c(object="QLearn"),
          definition = function(object, ...){ return(NULL) } )

setMethod(f = "FC",
          signature = c(object="QLearn"),
          definition = function(object, ...){ return(NULL) } )

setMethod(f = "FM",
          signature = c(object="OptimalSeq"),
          definition = function(object, ...){ return(NULL) } )

setMethod(f = "FC",
          signature = c(object="OptimalSeq"),
          definition = function(object, ...){ return(NULL) } )

setMethod(f = "FM",
          signature = c(object="OptimalClass"),
          definition = function(object, ...){ return(NULL) } )

setMethod(f = "FC",
          signature = c(object="OptimalClass"),
          definition = function(object, ...){ return(NULL) } )

setMethod(f = "fittedMain",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(FM(object)) } )

setMethod(f = "fittedCont",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(FC(object)) } )


setMethod(f = "genetic",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(Genetic(object)) } )

setMethod(f = "optTx",
          signature = c(x = "DynTxRegime", 
                        newdata = "missing"),
          definition = function(x, ...){ return(OptTx(x, ...)) } )

setMethod(f = "optTx",
          signature = c(x = "DynTxRegime", 
                        newdata = "data.frame"),
          definition = function(x, newdata, ...){ return(OptTx(x, newdata, ...)) } )

setMethod(f = "outcome",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(Outcome(object,...)) } )

setMethod(f = "plot",
          signature = c(x = "DynTxRegime"),
          definition = function(x, y, ...){ Plot(x, ...) } )

setMethod(f="print",
          signature=c(x="DynTxRegime"),
          definition = function(x){
                         Print(x)
                       } )
setMethod(f = "propen",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(Propen(object)) } )

setMethod(f = "DTRstep",
          signature=c(object="DynTxRegime"),
          definition = function(object,...){ return(Step(object)) } )

setMethod(f = "regimeCoef",
          signature = c(object="DynTxRegime"),
          definition = function (object, ...){ return(RegimeCoef(object)) } )

setMethod(f = "residuals",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(Residuals(object)) } )

setMethod(f="show",
          signature=c(object="DynTxRegime"),
          definition = function(object){
                         Show(object)
                       } )

setMethod(f = "stdDev",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(StdDev(object)) } )

setMethod(f = "summary",
          signature = c(object="DynTxRegime"),
          definition = function(object, ...){ return(MySummary(object, ...)) } )
