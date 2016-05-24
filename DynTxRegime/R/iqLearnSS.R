#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# iqLearnSS - Second stage in IQ-Learning algorithm                            #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain  : an object of class modelObj that defines the models and R methods  #
#           to be used to obtain parameter estimates and predictions for main  #
#           effects component of outcome regression.                           #
#           NULL is an acceptable value if moCont is defined.                  #
#                                                                              #
# moCont  : an object of class modelObj that defines the models and R methods  #
#           to be used to obtain parameter estimates and predictions for       #
#           contrasts component of outcome regression.                         #
#           NULL is an acceptable value if moMain is defined.                  #
#                                                                              #
# data    : data.frame of covariates and treatment histories                   #
#                                                                              #
# response: response vector                                                    #
#                                                                              #
# txName  : character string of treatment variable in data                     #
#                                                                              #
# iter    : an integer                                                         #
#                                                                              #
#           >=1 if moMain and moCont are to be fit iteratively                 #
#           The value is the maximum number of iterations.                     #
#           Note the iterative algorithm is as follows:                        #
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
#= Returns an object of class iqLearnSS.                                      =#
#=                                                                            =#
#==============================================================================#
iqLearnSS <- function(..., 
                      moMain,
                      moCont,
                      data,
                      response,
                      txName,
                      iter = 0L,
                      suppress = FALSE){

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Verify Input                         ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #--------------------------------------------------------------------------#
  # iter must be an integer                                                  #
  #--------------------------------------------------------------------------#
  if( !is(iter, "integer") ) iter <- as.integer(round(iter,0L))

  #--------------------------------------------------------------------------#
  # moMain must be either an object of class modelObj or NULL                #
  #--------------------------------------------------------------------------#
  if( !is(moMain, 'modelObj') && !is(moMain,'NULL') ) {
    UserError("input",
              "The class of 'moMain' must be one of {modelObj, NULL}.")
  }

  #--------------------------------------------------------------------------#
  # moCont must be either an object of class modelObj or NULL                #
  #--------------------------------------------------------------------------#
  if( !is(moCont, 'modelObj') && !is(moCont,'NULL') ) {
    UserError("input",
              "The class of 'moCont' must be one of {modelObj, NULL}.")
  }

  #--------------------------------------------------------------------------#
  # At least one of {moMain, moCont} must be an object of class modelObj     #
  # If either is NULL, iterative algorithm is not appropriate.               #
  #--------------------------------------------------------------------------#
  if( is(moMain, "NULL") && is(moCont, "NULL") ){
    UserError("input",
              "At least one of {moMain, moCont} must of provided.")
  } else if( is(moMain, "NULL") || is(moCont, "NULL") ) {
    if( iter != 0L ) {
      iter = 0L
      cat("Input variable iter reset to 0.\n")
    }
  }

  #--------------------------------------------------------------------------#
  # data must be an object of class data.frame                               #
  #--------------------------------------------------------------------------#
  if( !is(data, "data.frame") ) {
    UserError("input",
              "`data' must be a data.frame.")
  }

  #--------------------------------------------------------------------------#
  # response must be an object of class vector                               #
  #--------------------------------------------------------------------------#
  if( !is(response, "vector") ) {
    UserError("input",
              "`response' must be a vector.")
  }

  #--------------------------------------------------------------------------#
  # txName must be an object of class character                              #
  #--------------------------------------------------------------------------#
  if( !is(txName, "character") ) {
    UserError("input",
              "'txName' must be a character.")
  }

  #--------------------------------------------------------------------------#
  # Verify that the column exists                                            #
  #--------------------------------------------------------------------------#
  txVec <- try(data[,txName], silent = TRUE)

  if( is(txVec,"try-error") ) {
    UserError("input",
              paste(txName, " not found in data.", sep="") )
  }

  #--------------------------------------------------------------------------#
  # Treatment must be an integer.                                            #
  # If a factor, throw error. If numeric, recast as integer.                 #
  #--------------------------------------------------------------------------#
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

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Calculation                          ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #--------------------------------------------------------------------------#
  # Process treatment information.                                           #
  # Note that fSet cannot be specified for IQ-Learning methods and thus it   #
  # defaults to NULL in call to txProcess.                                   #
  #--------------------------------------------------------------------------#
  txInfo <- txProcess(txVar = txName, 
                      data = data, 
                      fSet = NULL)

  #--------------------------------------------------------------------------#
  # verify treatments are binary (-1L,1L)                                    #
  #--------------------------------------------------------------------------#
  sset <- SuperSet(txInfo)

  if( as.integer(round(length( sset ),0L)) != 2L ){
    UserError("input",
              "Only binary treatment options are allowed")
  }

  if( !all(sset %in% c("-1","1")) ) {
    UserError("input",
              "Treatment must be coded as {-1,1}")
  }

  #--------------------------------------------------------------------------#
  # Obtain fit of moMain and moCont models                                   #
  #--------------------------------------------------------------------------#
  est <- iqEst(moMain = moMain, 
               moCont = moCont, 
               data = data,
               response = response, 
               txInfo = txInfo, 
               iter = iter)

  #--------------------------------------------------------------------------#
  # Calculate Q-Functions at each tx                                         #
  #--------------------------------------------------------------------------#
  qFunctions <- matrix(data = 0.0,
                       nrow = nrow(data),
                       ncol = 2L,
                       dimnames = list(NULL,sset))

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
  # Create object of class iqLearnSS to return to user                       #
  #--------------------------------------------------------------------------#
  result <- new("IQLearnSS",
                txName = txName, 
                call = match.call(),
                qFunctions = qFunctions,
                est )

  if( !suppress ) show(result)

  return(result)
  
}
