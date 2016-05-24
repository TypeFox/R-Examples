#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                 CLASS IQEst                                  #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results for a q-learning step of the iq-learning algorithm          #
#                                                                              #
#  fitObj : a SimpleFit or IterateFit object.                                  #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQEst",
         slots = c( fitObj = "SimpleFit or IterateFit" ) )

setMethod(f = "Classif",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Base",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( Base(object@fitObj) ) 
                       } )

setMethod(f = "Coef", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( Coef(object = object@fitObj, ...) )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( FitObject(object = object@fitObj, ...) ) 
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( FittedCont(object = object@fitObj, ...) ) 
                       } )

setMethod(f = "FittedMain", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( FittedMain(object = object@fitObj, ...) ) 
                       } )

setMethod(f = "Genetic",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( ModelObjectFit(object = object@fitObj, ...) ) 
                       } )

setMethod(f = "MySummary", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){
                         return( MySummary(object = object@fitObj, ...) )
                       } )

setMethod(f = "Outcome",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         return( FitObject(object = object@fitObj, ...) ) 
                       } )

setMethod(f = "Plot", 
          signature = c(x = "IQEst"), 
          definition = function(x, suppress=FALSE, ...){ 
                         Plot(x = x@fitObj, suppress = suppress, ...) 
                       } )

setMethod(f = "PredictCont", 
          signature = c(object = "IQEst", 
                        newdata="data.frame"), 
          definition = function(object, newdata, ...){
                         res <- PredictCont(object = object@fitObj, 
                                            newdata = newdata, ...)
                         return( res )
                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "IQEst", 
                        newdata="data.frame"), 
          definition = function(object, newdata, ...){
                         res <- PredictMain(object = object@fitObj, 
                                            newdata = newdata, ...)
                         return( res )
                       } )

setMethod(f = "Print", 
          signature = c(x = "IQEst"), 
          definition = function(x, ...){ 
                         Print(x@fitObj, ...) 
                       } )

setMethod(f = "Propen",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "RegimeCoef",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Residuals",
          signature = c(object="IQEst"),
          definition = function(object, ...){ 
                         return( Residuals(object@fitObj, ...) ) 
                       } )

setMethod(f = "Show", 
          signature = c(object = "IQEst"), 
          definition = function(object, ...){ 
                         Show(object@fitObj, ...) 
                       } )

setMethod(f = "StdDev",    
          signature = c(object = "IQEst"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                  CLASS IQMin                                 #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQMin",
         slots = c(       call = "call",
                        txName = "character",
                    qFunctions = "matrix") )

setMethod(f = "Call",
          signature = c(object = "IQMin"),
          definition = function(object, ...){ 
                         return( object@call ) 
                       } )

setMethod(f = "TxName", 
          signature = c(object = "IQMin"),
          definition = function(object, ...){  
                         return( object@txName )  
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                CLASS IQLearnSS                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQLearnSS", 
         contains = c("IQEst", "IQMin", "DynTxRegime") )

setMethod(f = "Est",
          signature = c(x = "IQLearnSS"),
          definition = function(x, ...){
                         maxValue <- apply(x@qFunctions, 1L, max)
                         return( mean(maxValue) )
                       } )

setMethod(f = "OptTx",
          signature = c(x = "IQLearnSS", 
                        newdata = "missing"),
          definition = function(x, newdata, ...){

                         q2opt <- max.col(x@qFunctions, ties.method="first")
                         optTx <- as.integer(colnames(x@qFunctions)[q2opt])

                         return( list("qFunctions" = x@qFunctions,
                                       "optimalTx" = optTx) )
                       } )

setMethod(f = "OptTx",
          signature = c(x = "IQLearnSS", 
                        newdata = "data.frame"),
          definition = function(x, newdata, ...){

                         ps <- iqLearn_pm(object = x, 
                                          newdata = newdata)

                         colnames(ps) <- c("-1","1")

                         q2opt <- max.col(ps, ties.method="first")
                         optTx <- as.integer(colnames(ps)[q2opt])

                         return( list("qFunctions" = ps,
                                       "optimalTx" = optTx) )
                     } )


setMethod(f = "Print",
          signature = c(x = "IQLearnSS"),
          definition = function(x, ...){
            cat("\n")
            print(Call(x))
            cat("\n")
            Print(x@fitObj)
            cat("Mean of Value Function: ", Est(x), "\n\n")
          } )

setMethod(f = "Step", 
          signature = c(object = "IQLearnSS"), 
          definition = function(object, ...){ 
                         return("IQ-Learning: Second Stage") 
                       } )

setMethod(f = "Show",
          signature = c(object = "IQLearnSS"),
          definition = function(object, ...){
            Show(object@fitObj)
            cat("\nMean of Value Function: ", Est(object), "\n\n")
          } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                CLASS IQLearnFS                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Virtual class to combine all first stage results into a common class         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQLearnFS",
         contains = c("IQMin", "VIRTUAL", "DynTxRegime") )

setMethod(f = "Classif",    
          signature = c(object = "IQLearnFS"), 
          definition = function(object, ...){
                         return( NULL )
                       } )


setMethod(f = "Est",
          signature = c(x = "IQLearnFS"),
          definition = function(x, y, z, dens,...){

                         if( missing(dens) ) {
                           msg <- "dens must be one of {norm, nonpar}"
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( !(dens %in% c("norm", "nonpar")) ) {
                           msg <- "dens must be one of {norm, nonpar}"
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( missing(y) || missing(z) ) {
                           msg <- paste("Must provide objects for main effect,",
                                        "contrast, and variance steps.",sep="")
                           e <- simpleError(msg)
                           stop(e)
                         }

                         objs <- list(x,y,z)
                         classes <- c(class(x),class(y),class(z))

                         mo <- which(classes == "IQLearnFS_ME")
                         co <- which(classes == "IQLearnFS_C")
                         so <- which(classes == "IQLearnFS_VHet")

                         if( length(mo) == 0L ) {
                           msg <- "No main effects step found in input."
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( length(co) == 0L ) {
                           msg <- "No contrast step found in input."
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( length(so) == 0L ) {
                           so <- which(classes == "IQLearnFS_VHom")
                           if( length(so) == 0L ) {
                             msg <- "No variance step found in input."
                             e <- simpleError(msg)
                             stop(e)
                           }
                         }

                         x1 <- iqLearn_optTx1(mainObj = objs[[mo]], 
                                              cmObj = objs[[co]], 
                                              sigObj = objs[[so]], 
                                              dens = dens)

                         values <- apply(x1$qFunctions, 1L, max)

                         return( mean(values) )
                       } )

setMethod(f = "OptTx",
          signature = c(x = "IQLearnFS", 
                        newdata = "missing"),
          definition = function(x, newdata, y, z, dens,...){

                         if( missing(dens) ) {
                           msg <- "dens must be one of {norm, nonpar}"
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( !(dens %in% c("norm", "nonpar")) ) {
                           msg <- "dens must be one of {norm, nonpar}"
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( missing(y) || missing(z) ) {
                           msg <- paste("Must provide objects for main effect,",
                                        "contrast, and variance steps.",sep="")
                           e <- simpleError(msg)
                           stop(e)
                         }

                         objs <- list(x,y,z)
                         classes <- c(class(x),class(y),class(z))

                         mo <- which(classes == "IQLearnFS_ME")
                         co <- which(classes == "IQLearnFS_C")
                         so <- which(classes == "IQLearnFS_VHet")

                         if( length(mo) == 0L ) {
                           msg <- "No main effects step found in input."
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( length(co) == 0L ) {
                           msg <- "No contrast step found in input."
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( length(so) == 0L ) {
                           so <- which(classes == "IQLearnFS_VHom")
                           if( length(so) == 0L ) {
                             msg <- "No variance step found in input."
                             e <- simpleError(msg)
                             stop(e)
                           }
                         }

                         x1 <- iqLearn_optTx1(mainObj = objs[[mo]], 
                                              cmObj = objs[[co]], 
                                              sigObj = objs[[so]], 
                                              dens = dens)

                         return( x1 )
                       } )

setMethod(f = "OptTx",
          signature = c(x = "IQLearnFS", 
                        newdata = "data.frame"),
          definition = function(x, newdata, y, z, dens,...){

                         if( missing(dens) ) {
                           msg <- "dens must be one of {norm, nonpar}"
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( !(dens %in% c("norm", "nonpar")) ) {
                           msg <- "dens must be one of {norm, nonpar}"
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( missing(y) || missing(z) ) {
                           msg <- paste("Must provide objects for main effect,",
                                        "contrast, and variance steps.",sep="")
                           e <- simpleError(msg)
                           stop(e)
                         }

                         objs <- list(x,y,z)
                         classes <- c(class(x),class(y),class(z))

                         mo <- which(classes == "IQLearnFS_ME")
                         co <- which(classes == "IQLearnFS_C")
                         so <- which(classes == "IQLearnFS_VHet")

                         if( length(mo) == 0L ) {
                           msg <- "No main effects step found in input."
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( length(co) == 0L ) {
                           msg <- "No contrast step found in input."
                           e <- simpleError(msg)
                           stop(e)
                         }

                         if( length(so) == 0L ) {
                           so <- which(classes == "IQLearnFS_VHom")
                           if( length(so) == 0L ) {
                             msg <- "No variance step found in input."
                             e <- simpleError(msg)
                             stop(e)
                           }
                         }

                         x1 <- iqLearn_optTx1(mainObj = objs[[mo]], 
                                              cmObj = objs[[co]], 
                                              sigObj = objs[[so]], 
                                              dens = dens,
                                              newdata = newdata)

                         return( x1 )
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                              CLASS IQLearnFS_ME                              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by fitting models to second-stage main effect  in  #
# IQ-learning                                                                  #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQLearnFS_ME", 
         contains = c("IQEst", "IQLearnFS") )

setMethod(f = "Step", 
          signature = c(object = "IQLearnFS_ME"), 
          definition = function(object, ...){ 
                         return( "IQ-Learning: First Stage Regression of ME" ) 
                       } )

setMethod(f = "Print",
          signature = c(x = "IQLearnFS_ME"),
          definition = function(x, ...){
            cat("\n")
            print(Call(x))
            cat("\n")
            Print(x@fitObj)
          } )
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                              CLASS IQLearnFS_C                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by fitting models to second-stage contrast in      #
# IQ-learning                                                                  #
#                                                                              #
#   txVar     : vector of txs given at first-stage                             #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQLearnFS_C", 
         slots = c( txVec = "integer"),
         contains = c("IQEst", "IQLearnFS") )

setMethod(f = "Print",
          signature = c(x = "IQLearnFS_C"),
          definition = function(x, ...){
            cat("\n")
            print(Call(x))
            cat("\n")
            Print(x@fitObj)
          } )
setMethod(f = "TxVec",  
          signature = c(object="IQLearnFS_C"), 
          definition = function(object, ...){ 
                         return( object@txVec ) 
                       } )

setMethod(f = "Step", 
          signature = c(object = "IQLearnFS_C"), 
          definition = function(object, ...){ 
                         return( "IQ-Learning: First Stage Regression of C" ) 
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                            CLASS IQLearnFS_VHet                              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by fitting models to second-stage contrast in      #
# IQ-learning                                                                  #
#                                                                              #
#   scale : normalization for standard residual                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQLearnFS_VHet", 
         slots = c( scale = "numeric",
                    residuals = "numeric"),
         contains = c("IQEst", "IQLearnFS") )

setMethod(f = "qqPlot",
          signature = c(x="IQLearnFS_VHet"),
          definition = function(x, ...){
            x <- x@residuals
            qqnorm(x, ...)
            qqline(x)
          } )

setMethod(f = "Scale",  
          signature = c(object = "IQLearnFS_VHet"), 
          definition = function(object, ...){ 
                         return( object@scale ) 
                       } )
                 
setMethod(f = "Step", 
          signature = c(object = "IQLearnFS_VHet"), 
          definition = function(object, ...){ 
                         return( "IQ-Learning: First Stage Variance; Log-Linear" ) 
                       } )

setMethod(f = "Print",
          signature = c(x = "IQLearnFS_VHet"),
          definition = function(x, ...){
            cat("\n")
            print(Call(x))
            cat("\n")
            Print(x@fitObj)
          } )
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                            CLASS IQLearnFS_VHom                              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IQLearnFS_VHom", 
         slots = c(residuals = "numeric",
                      stdDev = "numeric"),
         contains = c("IQLearnFS") )

setMethod(f = "Coef",    
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "FitObject",    
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Genetic",    
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Outcome",    
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Print",
          signature = c(x = "IQLearnFS_VHom"),
          definition = function(x, ...){
            cat("\n")
            print(Call(x))
            cat("\n")
            cat("Standard Deviation: ", x@stdDev, "\n", sep="")
          } )

setMethod(f = "Propen",    
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "RegimeCoef",    
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Residuals", 
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){ 
                         return( object@residuals ) 
                       } )
  
setMethod(f = "Show",
          signature = c(object = "IQLearnFS_VHom"),
          definition = function(object, ...){
            cat("\n")
            show(Call(object))
            cat("\n")
            cat("Standard Deviation: ", object@stdDev, "\n", sep="")
          } )

setMethod(f = "StdDev", 
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){ 
                         return( object@stdDev ) 
                       } )
      
setMethod(f = "Step", 
          signature = c(object = "IQLearnFS_VHom"), 
          definition = function(object, ...){ 
                         return( "IQ-Learning: First Stage Variance; Constant" ) 
                       } )


