#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# iqLearnFSV - obtain the variance of the contrast mean                        #
#              This is the variance step of the iqLearn method.                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# For homoskedastic :                                                          #
#                                                                              #
# object  : object of class iqLearnFS_C                                        #
#                                                                              #
# For heteroskedastic :                                                        #
#                                                                              #
# object  : object of class iqLearnFS_C                                        #
#                                                                              #
# moMain  : an object of class modelObj that defines the models and R          #
#           methods to be used to obtain parameter estimates and predictions   #
#           for main effects component of outcome regression.                  #
#           NULL is an acceptable value if moCont is defined.                  #
#                                                                              #
# moCont  : an object of class modelObj that defines the models and R          #
#           methods to be used to obtain parameter estimates and predictions   #
#           for contrast component of outcome regression.                      #
#           NULL is an acceptable value if moMain is defined.                  #
#                                                                              #
# data    : data.frame of covariates and treatment histories                   #
#                                                                              #
# iter    : an integer                                                         #
#                                                                              #
#           >=1 if moMain and moCont are to be fitted iteratively              #
#           The value is the maximum number of iterations.                     #
#           Note the iterative algorithms is as follows:                       #
#           Y = Ymain + Ycont                                                  #
#            (1) hat(Ycont) = 0                                                #
#            (2) Ymain = Y - hat(Ycont)                                        #
#            (3) fit Ymain ~ moMain                                            #
#            (4) set Ycont = Y - hat(Ymain)                                    #
#            (5) fit Ycont ~ moCont                                            #
#            (6) Repeat steps (2) - (5) until convergence or                   #
#            a maximum of iter iterations.                                     #
#                                                                              #
#           <=0 moMain and moCont will be combined and fit as a single object. #
#                                                                              #
#           Either categorical or integer data can be provided for the tx.     #
#           If categorical, the fitted contrast and main effects are defined   #
#           relative to the base category {defined as levels()[1]}. The values #
#           may not be those returned by predict(object) if iterate fits are   #
#           used. If integer, the fitted contrast and main effects are defined #
#           relative to no tx (tx = 0).                                        #
#                                                                              #
#           Note that if iter <= 0, all non-model components of the            #
#           moMain and moCont must be identical                                #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class iqLearnFS_VHom or iqLearnFS_VHet                =#
#=                                                                            =#
#==============================================================================#
iqLearnFSV <- function(object, 
                       ...,
                       moMain = NULL, 
                       moCont = NULL, 
                       data = NULL, 
                       iter = 0L,
                       suppress = FALSE){

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Verify Input                         ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if( !is(object,'IQLearnFS_C') ){
    msg <- "'object' must be an object returned by a call to iqLearnFSC()."
    UserError("input", msg)
  }

  if( is.null(moMain) && is.null(moCont) ){

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    #++++++                       Calculation                        ++++++#
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    result <- iqLearnVarHom(object)

  } else {

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    #++++++                       Verify Input                       ++++++#
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    #----------------------------------------------------------------------#
    # iter must be an integer                                              #
    #----------------------------------------------------------------------#
    if( !is(iter, "integer") ) iter <- as.integer(round(iter,0L))

    #----------------------------------------------------------------------#
    # moMain must be either an object of class modelObj or NULL            #
    #----------------------------------------------------------------------#
    if( !is(moMain, 'modelObj') && !is(moMain,'NULL') ) {
      UserError("input",
                "The class of 'moMain' must be one of {modelObj, NULL}.")
    }

    #----------------------------------------------------------------------#
    # moCont must be either an object of class modelObj or NULL            #
    #----------------------------------------------------------------------#
    if( !is(moCont, 'modelObj') && !is(moCont,'NULL') ) {
      UserError("input",
                "The class of 'moCont' must be one of {modelObj, NULL}.")
    }

    #----------------------------------------------------------------------#
    # At least one of {moMain, moCont} must be an object of class modelObj #
    # If either is NULL, iterative algorithm is not appropriate.           #
    #----------------------------------------------------------------------#
    if( is(moMain,'NULL') && is(moCont,'NULL') ){
      UserError("input",
                "At least one of {moMain, moCont} must of provided.")
    } else if( is(moMain,'NULL') || is(moCont,'NULL') ) {
      if( iter != 0L ) {
        iter = 0L
        cat("Input variable iter reset to 0.\n")
      }
    }

    #----------------------------------------------------------------------#
    # data must be an object of class data.frame                           #
    #----------------------------------------------------------------------#
    if( !is(data, "data.frame") ) {
      UserError("input",
                "`data' must be a data.frame.")
    }

    #----------------------------------------------------------------------#
    # txName must be an object of class character                          #
    #----------------------------------------------------------------------#
    txName <- TxName(object)
    if( !is(txName, "character") ) {
      DeveloperError("'txName' must be a character.", "iqLearnVar - Het")
    }

    #----------------------------------------------------------------------#
    # Verify that the column exists                                        #
    #----------------------------------------------------------------------#
    txVec <- try(data[,txName], silent = TRUE)

    if( is(txVec,"try-error") ) {
      UserError("input",
                paste(txName, " not found in data.", sep="") )
    }

    #----------------------------------------------------------------------#
    # Treatment must be an integer.                                        #
    # If a factor, throw error. If numeric, recast as integer.             #
    #----------------------------------------------------------------------#
    if( is(txVec,"factor") ) {
        UserError("input",
                  "Treatment variable must be an integer.")
    } else {
      if( !isTRUE(all.equal(txVec, round(txVec,0L))) ) {
        UserError("input",
                  "Treatment variable must be an integer.")
      }
      data[,txName] <- as.integer(round(data[,txName],0L))
    }

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    #++++++                       Calculation                        ++++++#
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    result <- iqLearnVarHet(object = object, 
                            moMain = moMain,  
                            moCont = moCont,  
                            data = data,  
                            iter = iter)
  }

  if( !suppress ) show(result)
  return(result)
}


iqLearnVarHom <-  function(object){

  cmeanResids <- as.vector(Residuals(object))
  stdDev <- sd(cmeanResids)
  stdResids <- cmeanResids/stdDev

  result <- new("IQLearnFS_VHom",
                txName = TxName(object),
                call = match.call(),
                residuals = stdResids,
                stdDev = stdDev)

  return(result)
}


iqLearnVarHet <- function(object, 
                          moMain, 
                          moCont, 
                          data, 
                          iter){

  cmeanResids <- as.vector(Residuals(object))
  response <- log(cmeanResids*cmeanResids)

  txName <- TxName(object)

  #--------------------------------------------------------------------------#
  # Obtain fit of moMain and moCont models                                   #
  #--------------------------------------------------------------------------#
  est <- fitChoose(moMain = moMain, 
                   moCont = moCont, 
                   data = data,
                   response = response,
                   txName = txName, 
                   iter = iter)

  #--------------------------------------------------------------------------#
  # calculate the fitted response for patients included in fit of models     #
  #--------------------------------------------------------------------------#
  fitted <- FittedMain(est) + data[,txName]*FittedCont(est)

  #--------------------------------------------------------------------------#
  # Standardize the residuals of the contrast function after mean and        #
  # variance modeling                                                        #
  #--------------------------------------------------------------------------#
  sig <- exp(fitted/2.0)
  stdResids <- as.vector(cmeanResids/sig)
  sdr <- sd(stdResids)
  sd.stdResids <- 2.0*log(sdr)
  sig <- sig*(sdr)
  stdResids <- as.vector(cmeanResids)/sig

  #--------------------------------------------------------------------------#
  # Calculate Q-Functions at each tx                                         #
  #--------------------------------------------------------------------------#
  qFunctions <- matrix(data = 0.0,
                       nrow = nrow(data),
                       ncol = 2L,
                       dimnames = list(NULL,c("-1","1")))

  #--------------------------------------------------------------------------#
  # Set treatment to -1 for all samples.                                     #
  #--------------------------------------------------------------------------#
  data[,txName] <- -1L

  #--------------------------------------------------------------------------#
  # Obtain predicted main effects and contrast for this new dataset.         #
  # Note that the contrast model explicitly includes the treatment variable  #
  # and thus does not need to be included in the prediction expression.      #
  #--------------------------------------------------------------------------#
  qFunctions[,1L] <- PredictMain(object = est, newdata = data) +
                     PredictCont(object = est, newdata = data)

  #--------------------------------------------------------------------------#
  # Set treatment to 1 for all samples.                                      #
  #--------------------------------------------------------------------------#
  data[,txName] <- 1L

  #--------------------------------------------------------------------------#
  # Obtain predicted main effects and contrast for this new dataset.         #
  # Note that the contrast model explicitly includes the treatment variable  #
  # and thus does not need to be included in the prediction expression.      #
  #--------------------------------------------------------------------------#
  qFunctions[,2L] <- PredictMain(object = est, newdata = data) +
                     PredictCont(object = est, newdata = data)

  #--------------------------------------------------------------------------#
  # Create new iqLearnFS_VHet object to return to user                       #
  #--------------------------------------------------------------------------#
  result <- new("IQLearnFS_VHet",
                fitObj = est,
                txName = TxName(object),
                qFunctions = qFunctions,
                residuals = stdResids,
                call = match.call(),
                scale = sd.stdResids)

  return(result)
}


