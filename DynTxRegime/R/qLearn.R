#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# qLearn - public function to perform a step of the Q-Learning algorithm       #
#          If an object of class QLearn is passed, it is assumed to be the     #
#          preceding step of the Q-Learning algorithm and models are fit       #
#          using the Ytilde variable of the QLearn object. If a vector         #
#          is passed, it is assumed that this is the first step in the         #
#          Q-Learning algorithm and models are fit using the response.         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain  : an object of class modelObj or a list of objects of class          #
#           modelObjSubset, which define the models and R methods to           #
#           be used to obtain parameter estimates and predictions              #
#           for the main effects component of the outcome regression.          #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an acceptable value if moCont is defined.                  #
#                                                                              #
# moCont  : an object of class modelObj or a list of objects of class          #
#           modelObjSubset, which define the models and R methods to           #
#           be used to obtain parameter estimates and predictions              #
#           for the contrasts component of the outcome regression.             #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an acceptable value if moMain is defined.                  #
#                                                                              #
# data    : data frame of covariates and treatment histories                   #
#                                                                              #
# response: response vector or object of class qLearn from a previous          #
#           Q-Learning step.                                                   #
#                                                                              #
# txName  : character string giving column header of treatment variable        #
#           in data                                                            #
#                                                                              #
# fSet    : A function.                                                        #
#           This argument allows the user to specify the subset of tx          #
#           options available to a patient.                                    #
#           The functions should accept as input either                        #
#           1) explicit covariate names as given in column names of data       #
#           2) a vector of covariates (i.e. a row of a data.frame)             #
#           and must return a vector of tx options available to the            #
#           patient                                                            #
#           Note this function is used for an INDIVIDUAL patient, matrix       #
#           results are not appropriate                                        #
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
# base   : An integer indicating the base tx or NULL (ordinal tx)              #
#                                                                              #
# ...    : ignored                                                             #
#                                                                              #
# Function returns an object of class QLearn                                   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

qLearn <- function(...,
                   moMain,
                   moCont, 
                   data, 
                   response, 
                   txName, 
                   fSet = NULL, 
                   iter = 0L,
                   suppress = FALSE){


  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Verify Input                         ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #--------------------------------------------------------------------------#
  # data must be a data.frame object.                                        #
  #--------------------------------------------------------------------------#
  if( !is(data, 'data.frame') ) {

    UserError("input", 
              paste("'data' must be an object of class data.frame.", 
                    "Received an object of class ", 
                    paste(is(data), collapse=","), ".",
                    sep=""))

  }

  #--------------------------------------------------------------------------#
  # txName must be a character object.                                       #
  #--------------------------------------------------------------------------#
  if( !is(txName, 'character') ) {
    msg <- "'txName', must be a character object."
    UserError("input", msg)
  }

  #--------------------------------------------------------------------------#
  # Make sure treatment variable is in data.frame.                           #
  #--------------------------------------------------------------------------#
  txVec <- try(data[,txName], silent = TRUE)

  if( is(txVec,"try-error") ) {
    UserError("input",
              paste(txName, " not found in data.", sep="") )
  }

  #--------------------------------------------------------------------------#
  # Treatment must be either a factor or an integer.                         #
  #--------------------------------------------------------------------------#
  if( !is(txVec,"factor") ) {
    if( !isTRUE(all.equal(txVec, round(txVec,0L))) ) {
      UserError("input",
                "Treatment variable must be a factor or an integer.")
    }
    data[,txName] <- as.integer(round(data[,txName],0L))
  }

  #--------------------------------------------------------------------------#
  # Verify modeling objects.                                                 #
  #--------------------------------------------------------------------------#
  tempObj <- list()
  tempObj[[1]] <- moMain
  tempObj[[2]] <- moCont
  nms <- c("moMain","moCont")

  for (i in 1L:2L ) {

    if( is.list(tempObj[[i]]) ) {
      #------------------------------------------------------------------#
      # If obj is a list, verify verify that it contains information.    #
      #------------------------------------------------------------------#

      if( length(tempObj[[i]]) < 0.5 ) {
        UserError("input", 
                  paste(nms[i],"must have length > 0",sep="  "))
      }

      #------------------------------------------------------------------#
      # Verify that each element is a subset model object.               #
      #------------------------------------------------------------------#
      tst <- sapply(X = tempObj[[i]], FUN = class) != 'ModelObjSubset'
      if( any(tst) ) {
        UserError("input", 
                  paste("If class(",nms[i],") == list, ",
                        "all elements must be of class modelObjSubset.\n", 
                        "Received objects of class ", 
                        paste(is(tempObj[[i]]), collapse=", "), ".",
                        sep=""))
      }

      #------------------------------------------------------------------#
      # fSet must be provided when subset models are specified.          #
      #------------------------------------------------------------------#
      if( is(fSet, "NULL") ) {
        UserError("input", 
                  paste("When using objects of class modelObjSubset, ",
                        "fSet must be provided.", sep=""))
      }

      tempObj[[i]] <- new("ModelObjSubsetList",
                          loo = tempObj[[i]])

    } else if( !is(tempObj[[i]], "modelObj") && !is(tempObj[[i]], "NULL") ) {

      #------------------------------------------------------------------#
      # If single object is provided in obj, must be NULL or modelObj.   #
      #------------------------------------------------------------------#
      UserError("input", 
                paste("If modeling the superset of treatment options, ",
                       nms[i], "must be of class modelObj.\n", 
                      "Received an object of class ", 
                      paste(is(tempObj[[i]]),collapse=","), ".",
                      sep=""))

    }
  }

  moMain <- tempObj[[1L]]
  moCont <- tempObj[[2L]]


  #--------------------------------------------------------------------------#
  # At least one modeling object must be given.                              #
  #--------------------------------------------------------------------------#
  if( is(moMain, "NULL") && is(moCont, "NULL") ){
    UserError("input", 
              "Must provide at least one of {moMain, moCont}.")
  }

  #--------------------------------------------------------------------------#
  # Update Q-learning step information                                       #
  #--------------------------------------------------------------------------#
  if( is(response, "QLearn")  ){
    step <- IStep(response) + 1L
    response <- YTilde(response)
  } else if( is(response, "vector") ){
    step <- 1L
  } else {
    UserError("input", 
              paste("'response' must be a vector of responses or ",
                    "an object returned by a prior call to qLearn().",
                    sep=""))
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Calculation                          ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if( !suppress ) cat("\n\nStep ", step, " of Q-learning algorithm.\n\n")

  #--------------------------------------------------------------------------#
  # Process treatment option information.                                    #
  #--------------------------------------------------------------------------#
  if( is(fSet, "NULL") || is(fSet, "function") ){

    txInfo <- txProcess(txVar = txName, 
                        data = data, 
                        fSet = fSet)

  } else {

    UserError("input", 
              paste("If provided, fSet must be an object of class function.\n", 
                    "Received an object of class ", 
                    paste(is(fSet),collapse=","), ".",
                    sep=""))

  }

  #--------------------------------------------------------------------------#
  # Perform Q-Learning to obtain regression models.                          #
  #--------------------------------------------------------------------------#
  est <- qLearnEst(moMain = moMain, 
                   moCont = moCont,
                   data = data, 
                   response = response, 
                   iter = iter, 
                   txInfo = txInfo)

  #--------------------------------------------------------------------------#
  # Calculate Q-Functions at each tx                                         #
  #--------------------------------------------------------------------------#
  if( is(data[,txName],"factor") ) {
    nms <- levels(data[,txName])
  } else {
    nms <- SuperSet(txInfo)
  }

  qFunctions <- matrix(data = 0.0,
                       nrow = nrow(data),
                       ncol = length(nms),
                       dimnames = list(NULL, nms))

  for( i in 1L:length(nms) ) {

    if( is(data[,txName],"factor") ) {
      data[, txName] <- factor(rep(levels(data[,txName])[i],nrow(data)), 
                               levels = levels(data[,txName]))
    } else {
      data[, txName] <- as.integer(nms[i])
    }

    qFunctions[,i] <- PredictMain(object = est, newdata = data) +
                      PredictCont(object = est, newdata = data)

  }

  for( i in 1L:length(txInfo@subsets) ) {
    inss <- txInfo@ptsSubset %in% names(txInfo@subsets)[i]
    rmss <- !(nms %in% txInfo@subsets[[i]])
    qFunctions[inss,rmss] <- NA
  }

  fun <- function(x){
    if(all(is.na(x))) { 
      return(NA) 
    } else {
      return(which.max(x))
    }
  }

  q2opt <- apply(qFunctions,1,fun)

  if( is(data[,txName], "factor") ) {
    optTx <- factor(colnames(qFunctions)[q2opt],
                    levels = colnames(qFunctions))
  } else {
    optTx <- as.integer(colnames(qFunctions)[q2opt])
  }

  result <- new("QLearn", 
                step = step,
                qFunctions = qFunctions,
                optimalTx = optTx,
                call = match.call(),
                est,
                txInfo)

  if(!suppress) show(result)

  return(result)
}

