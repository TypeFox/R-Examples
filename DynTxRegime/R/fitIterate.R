#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# fitIterate - Fit main effects and contrasts models iteratively.              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain   : an object of class modelObj that defines the model and R methods  #
#            to be used to obtain parameter estimates and predictions for the  #
#            main effects component of an outcome regression step.             #
#            NULL is NOT an acceptable value for moMain.                       #
#                                                                              #
# moCont   : an object of class modelObj that defines the model and R methods  #
#            to be used to obtain parameter estimates and predictions for the  #
#            contrasts component of an outcome regression step.                #
#            NULL is NOT an acceptable value for moCont.                       #
#                                                                              #
# response : response vector                                                   #
#                                                                              #
# txName   : column header of data containing tx variable                      #
#                                                                              #
# data     : data.frame of covariates                                          #
#                                                                              #
# max.iter : maximum number of iterations                                      #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class iterateFit                                      =#
#=                                                                            =#
#==============================================================================#
fitIterate <- function(moMain,
                       moCont,
                       response,
                       txName,
                       data,
                       max.iter){

  tol <- 1.5e-8

  #--------------------------------------------------------------------------#
  # Subsets of data may be passed in. The following ensures that the         #
  # "levels" reflect the actual levels in the provided data.                 #
  #--------------------------------------------------------------------------#
  if( is(data[,txName], "factor") ) data[,txName] <- factor(data[,txName])

  #--------------------------------------------------------------------------#
  # convert contrast model to character string.                              #
  #--------------------------------------------------------------------------#
  tempCont <- terms(model(moCont))
  contPart <- paste(attr(tempCont,"term.labels"), collapse="+")

  #--------------------------------------------------------------------------#
  # Add treatment variable to contrast model.                                #
  #--------------------------------------------------------------------------#
  if( attr(tempCont,"intercept") > 0.5 ) {
    newForm <- paste("~ 0 + ", txName, " + ", 
                     txName, ":(", contPart, ")", sep="")
  } else if( attr(tempCont,"intercept") < 0.5 ){
    newForm <- paste("~ 0 + ", txName, ":(", contPart, ")", sep="")
  } 

  #--------------------------------------------------------------------------#
  # Convert character string to a formula object.                            #
  #--------------------------------------------------------------------------#
  newForm <- as.formula(newForm)

  #--------------------------------------------------------------------------#
  # Pull solver and predictor information from original moCont model obj     #
  #--------------------------------------------------------------------------#
  s <- solver(moCont)
  sa <- solverArgs(moCont)
  sa[[1L]] <- "formula"
  sa[[2L]] <- "data"
  p <- predictor(moCont)
  pa <- predictorArgs(moCont)
  pa[[1L]] <- "object"
  pa[[2L]] <- "newdata"

  #--------------------------------------------------------------------------#
  # Create new moCont model object using updated formula                     #
  #--------------------------------------------------------------------------#
  moCont <- buildModelObj(model = newForm, 
                          solver.method = s,
                          solver.args = sa,
                          predict.method = p,
                          predict.args = pa)

  #--------------------------------------------------------------------------#
  # Initialize contrast contribution to response to zero                     #
  #--------------------------------------------------------------------------#
  contrast <- matrix(data = 0.0, 
                     nrow = nrow(data), 
                     ncol = 1L)

  #--------------------------------------------------------------------------#
  # Iterate until the number of iterations reaches max.iter or difference    #
  # between iterations reaches tol                                           #
  #--------------------------------------------------------------------------#
  iter <- 0L
  while(TRUE){

    #----------------------------------------------------------------------#
    # Main effects component of response is the measured response less the #
    # current estimate for the contrast component.                         #
    #----------------------------------------------------------------------#
    YinternalY <- matrix(data = response-contrast,
                         ncol = 1L, 
                         dimnames = list(NULL,"YinternalY"))

    #----------------------------------------------------------------------#
    # Obtain fit for the main effects component                            #
    # (fit method from package modelObj)                                   #
    #----------------------------------------------------------------------#
    fit.main <- fit(object = moMain, 
                    data = data, 
                    response = YinternalY)

    #----------------------------------------------------------------------#
    # Calculate fitted main effects component of response                  #
    # (predict method from package modelObj)                               #
    #----------------------------------------------------------------------#
    main <- predict(object = fit.main, 
                    newdata = data)

    #----------------------------------------------------------------------#
    # Contrasts component of response is the measured response less the    #
    # current estimate for the main effects component.                     #
    #----------------------------------------------------------------------#
    YinternalY <- matrix(data = (response-main),
                         ncol = 1L, 
                         dimnames = list(NULL,"YinternalY"))

    #----------------------------------------------------------------------#
    # Obtain fit for the contrast component                                #
    # (fit method from package modelObj)                                   #
    #----------------------------------------------------------------------#
    fit.cont <- fit(object = moCont, 
                    data = data, 
                    response = YinternalY)

    #----------------------------------------------------------------------#
    # Calculate fitted contrasts component of the response                 #
    # (predict method from package modelObj)                               #
    #----------------------------------------------------------------------#
    contrast <- predict(object = fit.cont, 
                        newdata = data)

    #----------------------------------------------------------------------#
    # Compare the main effects and contrasts components of this iteration  #
    # to those of last iteration to determine max difference. If less      #
    # than tol, break loop.                                                #
    #----------------------------------------------------------------------#
    if( iter > 0L ){
      tst <- isTRUE(all.equal(contrast.old, contrast, tolerance=tol))
      tst <- isTRUE(all.equal(main.old, main, tolerance=tol)) && tst
      if( tst ) break
    }

    #----------------------------------------------------------------------#
    # Store current iteration values for comparison at next iteration      #
    #----------------------------------------------------------------------#
    contrast.old <- contrast
    main.old <- main
    iter <- iter + 1L
    #----------------------------------------------------------------------#
    # If hit max.iter without converfence, warn user, break from loop and  #
    # continue assuming convergence attained.                              #
    #----------------------------------------------------------------------#
    if( iter > max.iter ){
      warning("Convergence not attained within max iterations specified.")
      break
    }
  }

  #--------------------------------------------------------------------------#
  # Calculate the residual of the combined fit.                              #
  #--------------------------------------------------------------------------#
  residuals <- as.numeric(response - main - contrast)

  txVec <- data[,txName]

  if( is(data[,txName], "factor") ) {

    #----------------------------------------------------------------------#
    # If treatment is a factor variable:                                   #
    # Set tx to base case to obtain main effects contribution to predicted #
    # response. (predict method of package modelObj)                       #
    #----------------------------------------------------------------------#
    base <- levels(data[,txName])[1L]
    data[,txName] <- factor(rep(base,nrow(data)),
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
  mainBase <- predict(object = fit.main, 
                      newdata = data)

  #--------------------------------------------------------------------------#
  # Obtain fitted contrasts using 0 or base treatment.                       #
  #--------------------------------------------------------------------------#
  contrastBase <- predict(object = fit.cont, 
                          newdata = data)

  #--------------------------------------------------------------------------#
  # Calculate the fitted response for the base case.                         #
  #--------------------------------------------------------------------------#
  fittedBase <- mainBase + contrastBase

  #--------------------------------------------------------------------------#
  # The fitted contrast component is the full fitted response less the fitted#
  # response for the base case.                                              #
  # Remove the treatment variable if treatment is not a factor.              #
  #--------------------------------------------------------------------------#
  fittedCont <- (main + contrast) - fittedBase
  if( !is(txVec,"factor") ) {
    fittedCont <- fittedCont/txVec
    fittedCont[is.infinite(fittedCont)] <- 0.0
  }

  #--------------------------------------------------------------------------#
  # Create 'IterateFit' object to be returned to calling function.           #
  #--------------------------------------------------------------------------#
  result <- new("IterateFit",
                "baseLevel" = base,
                "txName" = txName,
                "modelObjectFitMain" = fit.main,
                "modelObjectFitCont" = fit.cont,
                "residuals" = residuals,
                "yMainHat" = as.numeric(fittedBase),
                "yContHat" = as.numeric(fittedCont))

  return(result)

}
