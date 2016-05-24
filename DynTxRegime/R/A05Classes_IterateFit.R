#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                              CLASS IterateFit                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by a call to modelObjFit with separate models      #
#                                                                              #
#   baseLevel          : character or numeric value of base treatment level    #
#                                                                              #
#   txName             : column header of data containing treatment variable   #
#                                                                              #
#   modelObjectFitMain : An object of class modelObjFit for main effects model #
#                                                                              #
#   modelObjectFitCont : An object of class modelObjFit for contrast model     #
#                                                                              #
#   residuals          : residuals of the combined fit                         #
#                                                                              #
#   yMainHat           : fitted main effects                                   #
#                                                                              #
#   yContHat           : fitted contrast                                       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Available Methods:                                                           #
#                                                                              #
# Base        : Returns value in slot `baseLevel'                              #
#                                                                              #
# Coef        : Returns the coefficients, as defined by the regression method. #
#               A list is returned with element names "MainEffect" & "Contrast"#
#                                                                              #
# FitObject   : Returns the value object returned by the regression method.    #
#               A list is returned with element names "MainEffect" & "Contrast"#
#                                                                              #
# FittedCont  : Returns value in slot 'yContHat'                               #
#                                                                              #
# FittedMain  : Returns value in slot 'yMainHat'                               #
#                                                                              #
# FitType     : Returns value in slot 'fitType'                                #
#                                                                              #
# ModelObjectFit : Returns value in slot 'modelObjectFit'                      #
#               A list is returned with element names "MainEffect" & "Contrast"#
#                                                                              #
# MySummary   : Returns a list of summary objects as defined by the regression #
#               method.                                                        #
#               A list is returned with element names "MainEffect" & "Contrast"#
#                                                                              #
# Plot        : Generates plots with augmented titles indicating the modeling  #
#               objects passed in  (MainEffect, Contrast)                      #
#                                                                              #
# PredictCont : Takes new data and uses predict function corresponding to      #
#               the regression method to obtain predictions for contrast       #
#               component. Uses base level.                                    #
#                                                                              #
# PredictMain : Takes new data and uses predict function corresponding to      #
#               the regression method to obtain predictions for main effect    #
#               component. Uses base level.                                    #
#                                                                              #
# Print       : Prints modelObjectFit                                          #
#                                                                              #
# Residuals   : Returns value in slot 'residuals'                              #
#                                                                              #
# Show        : Print modelObjectFit                                           #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IterateFit", 
         slots = c(baseLevel = "character or numeric",
                   txName = "character",
                   modelObjectFitMain = "modelObjFit",
                   modelObjectFitCont = "modelObjFit",
                   residuals = "numeric",
                   yMainHat = "numeric", 
                   yContHat = "numeric"))

setClass(Class = "IterateFitList",
         contains = "dpList")

setMethod(f = "Base",    
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ 
                         return( object@baseLevel ) 
                       } )

setMethod(f = "Coef",
          signature = c(object = "IterateFit"),
          definition = function(object, ...) {
                         me <- coef(object@modelObjectFitMain, ...)
                         cn <- coef(object@modelObjectFitCont, ...)
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){
                         me <- fitObject(object@modelObjectFitMain, ...)
                         cn <- fitObject(object@modelObjectFitCont, ...)
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ return( object@yContHat ) } )

setMethod(f = "FittedMain", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ return( object@yMainHat ) } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){
                         me <- object@modelObjectFitMain
                         cn <- object@modelObjectFitCont
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "MySummary",
          signature = c(object = "IterateFit"),
          definition = function(object, ...) {
                         me <- MySummary(object@modelObjectFitMain, ...)
                         cn <- MySummary(object@modelObjectFitCont, ...)
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "Plot", 
          signature = c(x = "IterateFit"), 
          definition = function(x, suppress=FALSE, ...){

                         argList <- list(...)
                         if( !suppress ) {
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- "MainEffect"
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- "MainEffect"
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (MainEffect)", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@modelObjectFitMain
                         do.call(what = plot, args = argList)

                         argList <- list(...)
                         if( !suppress ) {
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- "Contrast"
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- "Contrast"
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (Contrast)", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@modelObjectFitCont
                         do.call(what = plot, args = argList)
                       } )


setMethod(f = "PredictCont", 
          signature = c(object = "IterateFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         if( is(object@baseLevel, "character") ) {

                           newdata[,object@txName] <- 
                             as.factor(newdata[,object@txName])

                         } else {
                           newdata[,object@txName] <- as.integer(newdata[,object@txName])
                         }

                         main <- predict(object = object@modelObjectFitMain,  
                                         newdata = newdata)

                         contrast <- predict(object = object@modelObjectFitCont, 
                                             newdata = newdata)

                         fittedY <- main + contrast

                         n <- nrow(newdata)

                         if( is(object@baseLevel, "character") ) {

                           newdata[,object@txName] <- 
                             as.factor(rep(object@baseLevel,n))

                         } else {
                           newdata[,object@txName] <- object@baseLevel
                         }

                         mainBase <- predict(object = object@modelObjectFitMain,  
                                             newdata = newdata)

                         contrastBase <- predict(object = object@modelObjectFitCont, 
                                                 newdata = newdata)

                         fittedBase <- mainBase + contrastBase

                         fittedCont <- fittedY - fittedBase

                         return( fittedCont )
                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "IterateFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         n <- nrow(newdata)

                         if( is(object@baseLevel, "character") ) {
                           newdata[,object@txName] <- 
                             as.factor(rep(object@baseLevel,n))
                         } else {
                           newdata[,object@txName] <- object@baseLevel
                         }

                         mainBase <- predict(object = object@modelObjectFitMain,  
                                             newdata = newdata)

                         contrastBase <- predict(object = object@modelObjectFitCont, 
                                                 newdata = newdata)

                         fittedBase <- mainBase + contrastBase

                         return( fittedBase )
                       } )

setMethod(f = "Print", 
          signature = c(x = "IterateFit"), 
          definition = function(x, ...){
                         cat("\n *** MainEffect Fit *** \n")
                         show(x@modelObjectFitMain)
                         cat("\n *** Contrast Fit *** \n")
                         show(x@modelObjectFitCont)
                       } )

setMethod(f = "Residuals", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ 
                         return( object@residuals ) 
                       } )

setMethod(f = "Show", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){
                         cat("\n *** MainEffect Fit *** \n")
                         show(object@modelObjectFitMain)
                         cat("\n *** Contrast Fit *** \n")
                         show(object@modelObjectFitCont)
                       } )

if(!isClass("SimpleFit or IterateFit")){
  setClassUnion("SimpleFit or IterateFit", 
                members = c("SimpleFit","IterateFit"))
}


