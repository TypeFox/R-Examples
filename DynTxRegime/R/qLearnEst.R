#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# qLearnEst - internal function to perform a step of the Q-Learning algorithm  #
#             it is assumed that all input is appropriate for the call.        #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# moMain  : an object of class modelObj or modelObjSubsetList, which define    #
#           the models and R methods to be used to obtain parameter estimates  #
#           and predictions for the main effects component of the outcome      #
#           regression.                                                        #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an acceptable value if moCont is defined.                  #
#                                                                              #
# moCont  : an object of class modelObj or modelObjSubsetList, which define    #
#           the models and R methods to be used to obtain parameter estimates  #
#           and predictions for the contrast component of the outcome          #
#           regression.                                                        #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an acceptable value if moMain is defined.                  #
#                                                                              #
# data    : data frame of covariates and treatment histories                   #
#                                                                              #
# response: response vector                                                    #
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
#           Note that if iter <= 0, all non-model components of the            #
#           moMain and moCont must be identical                                #
#                                                                              #
# txInfo   : object of class TxInfo                                            #
#                                                                              #
# base     : An integer indicating the base tx or NULL (ordinal tx)            #
#                                                                              #
# ...    : ignored                                                             #
#                                                                              #
# Returns an object of class QLearnEst                                         #
################################################################################
qLearnEst <- function(moMain,
                      moCont,
                      data, 
                      response, 
                      txInfo,
                      iter, ...){

  optTx <- numeric(length=nrow(data))

  Ytilde <- as.numeric(response)

  #--------------------------------------------------------------------------#
  # Extract subset groups and vector indicating group assigned to each pt    #
  #--------------------------------------------------------------------------#
  subsets <- Subsets(txInfo)
  ptsSubset <- PtsSubset(txInfo)

  #--------------------------------------------------------------------------#
  # Pick a model object for purposes of identifying class and preprocessing  #
  #--------------------------------------------------------------------------#
  if( !is(moMain, "NULL") ) {
    tempObj <- moMain
  } else if( !is(moCont, "NULL") ) {
    tempObj <- moCont
  } else {
    DeveloperError("moMain and moCont are NULL.",
                   "qLearnEst")
  }

  if( is(tempObj,'ModelObjSubsetList' ) ) {

    dp <- DecisionPoint( tempObj[[1]] )

    fitObj <- list()

    for( j in 1L:length(tempObj) ){
      #------------------------------------------------------------------#
      # Retrieve nickname of subset for this model                       #
      #------------------------------------------------------------------#
      modelSubset <- Subset( tempObj[[j]] )

      #------------------------------------------------------------------#
      # Match nickname to list of subsets defined by fSet                #
      #------------------------------------------------------------------#
      isubset <- which( names(subsets) %in% modelSubset )

      if( length(isubset) < 0.5 ) {
        UserError("input",
                  paste("Unable to match subset ", 
                        paste(modelSubset,collapse=", "), 
                        " to a subset defined by fSet.", sep=""))
      }

      #------------------------------------------------------------------#
      # Use only patients whose treatment options match this group       #
      #------------------------------------------------------------------#
      use4fit <- ptsSubset %in% names(subsets)[isubset]

      if( sum(use4fit) < 0.5 ) {
        UserError("input",
                  paste("No observations match ", modelSubset, " at dp ", 
                  dp, ".", sep = ""))
      }

      #------------------------------------------------------------------#
      # Identify all required model objects                              #
      #------------------------------------------------------------------#
      if( !is(moMain, "NULL") && is(moCont, "NULL") ) {
        moMainTemp <- ModelObject(moMain[[j]])
        moContTemp <- NULL
      } else if( is(moMain, "NULL") && !is(moCont, "NULL") ) {
        moMainTemp <- NULL
        moContTemp <- ModelObject(moCont[[j]])
      } else if( !is(moMain, "NULL") && !is(moCont, "NULL") ) {
        moMainTemp <- ModelObject(moMain[[j]])
        k <- 1L
        found <- FALSE
        while( k <= length(moCont) ) {

          txgC <- Subset( moCont[[k]] )

          if( all(txgC %in% modelSubset) && all(modelSubset %in% txgC) ) {
            moContTemp <- ModelObject(moCont[[k]])
            found <- TRUE
            break
          }

          k <- k + 1L
        }
        if( !found ){
          UserError("input",
                    paste("Could not find matching moCont model for ", 
                          modelSubset, " at dp ", dp, ".", sep = ""))
        }
      }

      #------------------------------------------------------------------#
      # Fit model                                                        #
      #------------------------------------------------------------------#
      fitO <- fitChoose(moMain = moMainTemp, 
                        moCont = moContTemp, 
                        data = data[use4fit,,drop=FALSE], 
                        response = response[use4fit], 
                        txName = TxName(txInfo),
                        iter = iter)

      #------------------------------------------------------------------#
      # Calculate Q-functions                                            #
      #------------------------------------------------------------------#
      for( i in 1L:length(isubset) ){
        options <- subsets[[i]]

        #--------------------------------------------------------------#
        # Use only patients whose treatment options match this group   #
        #--------------------------------------------------------------#
        use4Q <- ptsSubset %in% names(subsets)[isubset[i]]

        vals <- matrix(data = 0.0, 
                       nrow = sum(use4Q), 
                       ncol = length(options),
                       dimnames = list(NULL,options))

        for( k in 1L:length(options) ){
          newdata <- data[use4Q,,drop=FALSE]
          if( is(data[,TxName(txInfo)], "factor") ) {
            newdata[,TxName(txInfo)] <- factor(options[k],
                                               levels = subsets[[i]])
          } else {
            newdata[,TxName(txInfo)] <- as.integer(options[k])
          }
          vals[,k] <- PredictMain(fitO, newdata) + PredictCont(fitO, newdata)
        }

        #----------------------------------------------------------------#
        # Determine which of the txs leads to the maximum value          #
        #----------------------------------------------------------------#
        cols <- max.col(vals, ties.method = "first")
        optTx[use4Q] <- options[cols[use4Q]]

        #----------------------------------------------------------------#
        # Store value function                                           #
        #----------------------------------------------------------------#
        Ytilde[use4Q] <- apply(X = vals,
                               MARGIN = 1L,
                               FUN = max)
      }

      fitObj[[j]] <- new("SubsetFit",
                         subset = modelSubset,
                         fitObj = fitO)

    }

    fitO <- new("SubsetFitList",
                decisionPoint = dp,
                txInfo = txInfo,
                loo = fitObj)

  } else if( is(tempObj,'modelObj') ) {

    #----------------------------------------------------------------------#
    # Eliminate patients with only 1 tx option from dataset for fit        #
    #----------------------------------------------------------------------#
    use4fit <- eliminateSingleTx(subsets=subsets, ptsSubset=ptsSubset)

    modelSubset <- SuperSet(txInfo)

    if( !is.null(moMain) && is.null(moCont) ) {
      moMainTemp <- moMain
      moContTemp <- NULL
    } else if( is.null(moMain) && !is.null(moCont) ) {
      moMainTemp <- NULL
      moContTemp <- moCont
    } else if( !is.null(moMain) && !is.null(moCont) ) {
      moMainTemp <- moMain
      moContTemp <- moCont
    }

    fitO <- fitChoose(moMain = moMainTemp, 
                      moCont = moContTemp, 
                      data = data[use4fit,,drop=FALSE], 
                      response = response[use4fit], 
                      txName = TxName(txInfo), 
                      iter = iter)

    #----------------------------------------------------------------------#
    # Use fit to calculate predicted response for each tx option           #
    #----------------------------------------------------------------------#
    vals <- matrix(data = 0.0, 
                   nrow = sum(use4fit), 
                   ncol = length(modelSubset),
                   dimnames = list(NULL,modelSubset))

    for( i in 1L:length(modelSubset) ){
      newdata <- data[use4fit,,drop=FALSE]
      if( is(data[,TxName(txInfo)], "factor") ) {
        newdata[,TxName(txInfo)] <- factor(modelSubset[i],
                                           levels = levels(data[,TxName(txInfo)]))
      } else {
        newdata[,TxName(txInfo)] <- as.integer(modelSubset[i])
      }
      vals[,i] <- PredictMain(fitO, newdata) + PredictCont(fitO, newdata)
    }

    for( i in 1L:length(subsets) ) {
      inss <- ptsSubset[use4fit] %in% names(subsets)[i]
      rmss <- !(modelSubset %in% subsets[[i]])
      vals[inss,rmss] <- NA
    }

    #----------------------------------------------------------------------#
    # Store value function                                                 #
    #----------------------------------------------------------------------#
    Ytilde[use4fit] <- apply(X = vals, 
                             MARGIN = 1L, 
                             FUN = max,
                             na.rm = TRUE)
  }

  #--------------------------------------------------------------------------#
  # calculate the residuals                                                  #
  # remove residuals for patients not used to obtain fits                    #
  #--------------------------------------------------------------------------#
  residuals <- response - Ytilde
  residuals <- residuals[!use4fit]

  result <- new("QLearnEst",
                fitObj = fitO,
                valueFunc = Ytilde)

  return(result)
}


