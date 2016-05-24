#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# fitCombined - Combines main effects and contrasts models to obtain           #
#               parameter estimates in a single call to fit method.            #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain   : an object of class modelObj that defines the model and R methods  #
#            to be used to obtain parameter estimates and predictions for the  #
#            main effects component of an outcome regression step.             #
#            NULL is an acceptable value for moMain if moCont is defined.      #
#            However, at least one of {moMain,moCont} must be defined.         #
#                                                                              #
# moCont   : an object of class modelObj that defines the model and R methods  #
#            to be used to obtain parameter estimates and predictions for the  #
#            contrasts component of an outcome regression step.                #
#            NULL is an acceptable value for moMain if moCont is defined.      #
#            However, at least one of {moMain,moCont} must be defined.         #
#                                                                              #
# response : response vector                                                   #
#                                                                              #
# txName   : column header of data containing tx variable                      #
#                                                                              #
# data     : data.frame of covariates                                          #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class SimpleFit                                       =#
#=                                                                            =#
#==============================================================================#
fitCombined <- function(moMain,
                        moCont,
                        response,
                        txName,
                        data){

  #--------------------------------------------------------------------------#
  # Recast response as matrix with column header YinternalY                  #
  #--------------------------------------------------------------------------#
  YinternalY <- matrix(data = response, 
                       ncol = 1L,  
                       dimnames = list(NULL,"YinternalY"))

  #--------------------------------------------------------------------------#
  # Subsets of data may be passed in. The following ensures that the         #
  # "levels" reflect the actual levels in the provided data.                 #
  #--------------------------------------------------------------------------#
  if( is(data[,txName], "factor") ) data[,txName] <- factor(data[,txName])

  #--------------------------------------------------------------------------#
  # Combine main effects and contrast formula objects                        #
  #--------------------------------------------------------------------------#
  if( !is(moCont, "NULL") && !is(moMain, "NULL") ){

    #----------------------------------------------------------------------#
    # If both moMain and moCont are provided, combine into a single        #
    # formula and create new modeling object.                              #
    #----------------------------------------------------------------------#
    type <- "Combined"

    #----------------------------------------------------------------------#
    # convert contrast model to character string.                          #
    #----------------------------------------------------------------------#
    tempCont <- terms(model(moCont))
    contPart <- paste(attr(tempCont,"term.labels"), collapse="+")

    #----------------------------------------------------------------------#
    # convert main effect model to character string.                       #
    #----------------------------------------------------------------------#
    tempMain <- terms(model(moMain))
    mainPart <- paste(attr(tempMain,"term.labels"), collapse="+")

    #----------------------------------------------------------------------#
    # If no intercept is in main effects model, insert 0                   #
    #----------------------------------------------------------------------#
    if( attr(tempMain,"intercept") < 0.5 ) {
      mainPart <- paste("0 + ", mainPart, sep="")
    }

    #----------------------------------------------------------------------#
    # Combine character strings to create new formula. Explicitly handle   #
    # intercept specification of contrast function, which introduces a     #
    # treatment term in the combined model.                                #
    #----------------------------------------------------------------------#
    if( attr(tempCont,"intercept") > 0.5 ) {
      newForm <- paste("~", mainPart, "+ ", txName, " + ",
                       txName, ":(", contPart, ")", sep="")
    } else if( attr(tempCont,"intercept") < 0.5 ) {
      newForm <- paste("~", mainPart, "+ ",
                       txName, ":(", contPart, ")", sep="")
    } 

    #----------------------------------------------------------------------#
    # Convert character string to a formula object.                        #
    #----------------------------------------------------------------------#
    newForm <- as.formula(newForm)

    #----------------------------------------------------------------------#
    # Pull solver and predictor information from original moMain model obj #
    #----------------------------------------------------------------------#
    s <- solver(moMain)
    sa <- solverArgs(moMain)
    sa[[1L]] <- "formula"
    sa[[2L]] <- "data"
    p <- predictor(moMain)
    pa <- predictorArgs(moMain)
    pa[[1L]] <- "object"
    pa[[2L]] <- "newdata"

    #----------------------------------------------------------------------#
    # Create new model object using updated formula                        #
    #----------------------------------------------------------------------#
    obj <- buildModelObj(model = newForm, 
                         solver.method = s,
                         solver.args = sa,
                         predict.method = p,
                         predict.args = pa)

  } else if( !is(moCont, "NULL") && is(moMain, "NULL") ) {

    #----------------------------------------------------------------------#
    # If only a contrast function is provided, add treatment variable      #
    #----------------------------------------------------------------------#
    type <- "Contrast"

    #----------------------------------------------------------------------#
    # convert contrast model to character string.                          #
    #----------------------------------------------------------------------#
    tempCont <- terms(model(moCont))
    contPart <- paste(attr(tempCont,"term.labels"), collapse="+")

    #----------------------------------------------------------------------#
    # Combine character strings to create new formula. Explicitly handle   #
    # intercept specification of contrast function, which introduces a     #
    # treatment term.                                                      #
    #----------------------------------------------------------------------#
    if( attr(tempCont,"intercept") > 0.5 ) {
      newForm <- paste("~ 0 + ", txName, " + ", 
                       txName, ":(", contPart, ")", sep="")
    } else if( attr(tempCont,"intercept") < 0.5 ){
      newForm <- paste("~ 0 + ", txName, ":(", contPart, ")", sep="")
    } 

    #----------------------------------------------------------------------#
    # Convert character string to a formula object.                        #
    #----------------------------------------------------------------------#
    newForm <- as.formula(newForm)

    #----------------------------------------------------------------------#
    # Pull solver and predictor information from original moCont model obj #
    #----------------------------------------------------------------------#
    s <- solver(moCont)
    sa <- solverArgs(moCont)
    sa[[1L]] <- "formula"
    sa[[2L]] <- "data"
    p <- predictor(moCont)
    pa <- predictorArgs(moCont)
    pa[[1L]] <- "object"
    pa[[2L]] <- "newdata"

    #----------------------------------------------------------------------#
    # Create new model object using updated formula                        #
    #----------------------------------------------------------------------#
    obj <- buildModelObj(model = newForm, 
                         solver.method = s,
                         solver.args = sa,
                         predict.method = p,
                         predict.args = pa)

  } else if( is(moCont, "NULL") && !is(moMain, "NULL")){

    #----------------------------------------------------------------------#
    # If only a main effects function is provided, do nothing.             #
    #----------------------------------------------------------------------#
    type <- "MainEffect"

    obj <- moMain

  } else {

    DeveloperError("moMain and moCont made it in as null.", "fitCombined")

  }

  #--------------------------------------------------------------------------#
  # Obtain fit                                                               #
  # (fit method of package modelObj)                                         #
  #--------------------------------------------------------------------------#
  fitCombined <- fit(object = obj, 
                     data = data, 
                     response = YinternalY)

  #--------------------------------------------------------------------------#
  # Obtain fitted response.                                                  #
  #--------------------------------------------------------------------------#
  fittedY <- predict(object = fitCombined,
                     newdata = data)

  n <- nrow(data)
  txVec <- data[,txName]

  if( is(data[,txName], "factor") ) {

    #----------------------------------------------------------------------#
    # If treatment is a factor variable:                                   #
    # Set tx to base case to obtain main effects contribution to predicted #
    # response. (predict method of package modelObj)                       #
    #----------------------------------------------------------------------#
    base <- levels(data[,txName])[1L]
    data[,txName] <- factor(rep(base,n),
                            levels = levels(data[,txName]))
  } else {

    #----------------------------------------------------------------------#
    # If treatment is numeric/integer                                      #
    # Set tx to zero  to obtain main effects contribution to predicted     #
    # response. (predict method of package modelObj)                       #
    #----------------------------------------------------------------------#
    base <- 0L
    data[,txName] <- base

  }

  #--------------------------------------------------------------------------#
  # Obtain fitted main effects using 0 or base treatment.                    #
  #--------------------------------------------------------------------------#
  fittedMain <- predict(object = fitCombined, 
                        newdata = data)

  #--------------------------------------------------------------------------#
  # Calculate fitted contrast, removing the treatment variable if treatment  #
  # is not a factor.                                                         #
  #--------------------------------------------------------------------------#
  fittedCont <- (fittedY - fittedMain)
  if( !is(txVec,"factor") ) {
    fittedCont <- fittedCont/txVec
    fittedCont[is.infinite(fittedCont)] <- 0.0
  }

  #--------------------------------------------------------------------------#
  # Calculate residuals.                                                     #
  #--------------------------------------------------------------------------#
  residuals <- as.numeric(YinternalY - fittedY)

  #--------------------------------------------------------------------------#
  # Create 'SimpleFit' object to be returned to calling function.            #
  #--------------------------------------------------------------------------#
  result <- new("SimpleFit",
                "baseLevel" = base,
                "txName" = txName,
                "fitType" = type, 
                "modelObjectFit" = fitCombined,
                "residuals" = residuals,
                "yMainHat" = as.vector(fittedMain),
                "yContHat" = as.vector(fittedCont))

  return(result)
}
