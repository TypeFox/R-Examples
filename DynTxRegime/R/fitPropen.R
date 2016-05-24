#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# fitPropen : Fit propensity score models.                                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moPropen      : an object of class modelObj or modelObjSubsetList, which     #
#                 define the models and R method to be used to obtain          #
#                 parameter estimates and predictions for the propensity of    #
#                 treatment.                                                   #
#                                                                              #
# txInfo        : an object of class txInfo                                    #
#                                                                              #
# data          : full data frame of covariates.                               #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class propenFit                                       =#
#=                                                                            =#
#==============================================================================#
fitPropen <- function(moPropen,
                      txInfo,
                      data){

  #--------------------------------------------------------------------------#
  # Retrieve feasible set information from treatment object.                 #
  #--------------------------------------------------------------------------#
  subsets <- Subsets(txInfo)
  ptsSubset <- PtsSubset(txInfo)
  txName <- TxName(txInfo)

  if( is(moPropen, "ModelObjSubsetList") ){
    #----------------------------------------------------------------------#
    # If subset modeling is used, must match models to subsets.            #
    #----------------------------------------------------------------------#

    propenFitObj <- list()

    for( k in 1L:length(moPropen) ){

      #------------------------------------------------------------------#
      # Extract the subset name(s) for the kth model                     #
      #------------------------------------------------------------------#
      modelSubset <- Subset(moPropen[[k]])

      #------------------------------------------------------------------#
      # Identify which fSet subsets match the model subsets              #
      # Note that a model can be applied to multiple subsets.            #
      #------------------------------------------------------------------#
      indSubset <- which( names(subsets) %in% modelSubset )

      if( length(indSubset) < 0.5 ) {
        UserError("input", 
                  paste("Unable to match moPropen subset ", modelSubset, 
                        " to a subset defined by the fSet.", sep=""))
      }

      #------------------------------------------------------------------#
      # Determine which samples are associated with the model subsets    #
      #------------------------------------------------------------------#
      use4fit <- ptsSubset %in% modelSubset

      if( sum(use4fit) < 0.5 ) {
        UserError("input", 
                  paste("No observations match moPropen model subset:",
                        paste(modelSubset,collapse=","),sep=" "))
      }

      #------------------------------------------------------------------#
      # Determine if the smallest or largest treatment is not present    #
      # in returned predictions. Note that some methods return           #
      # predictions for all treatment options.                           #
      #------------------------------------------------------------------#
      predArgs <- predictorArgs(moPropen[[k]]@modelObject)

      if( is(predArgs[["propen.missing"]], "NULL") ) {
        small <- TRUE
      } else if( tolower(predArgs[["propen.missing"]]) == "smallest" ) {
        small <- TRUE
        predArgs[["propen.missing"]] <- NULL
      } else if( tolower(predArgs[["propen.missing"]]) == "largest" ) {
        small <- FALSE
        predArgs[["propen.missing"]] <- NULL
      } 
      predictorArgs(moPropen[[k]]@modelObject) <- predArgs
      
      #------------------------------------------------------------------#
      # Subset modeling causes problems with identifiability of tx levels#
      # The following tries to track what is being estimated. If tx is   #
      # a factor, the tx vector for the subset of data is recast as      #
      # another factor variable to get the correct levels. If tx is an   #
      # integer, the unique tx values for the subset of data are sorted  #
      # It is thus assumed that the prediction method will return the    #
      # probability matrix in the default order of factors or in sorted  #
      # order of the integer values.                                     #
      #------------------------------------------------------------------#
      tData <- data[use4fit,,drop=FALSE]

      if( is.factor(data[,txName]) ) {
        #--------------------------------------------------------------#
        # If treatment is a factor, pull cases that are to be fit and  #
        # reset levels                                                 #
        #--------------------------------------------------------------#
        tempTx <- factor(data[use4fit,txName])
        levs <- levels(tempTx)
      } else {
        #--------------------------------------------------------------#
        # If treatment is not a factor, pull cases that are to be fit  #
        # and define levels as the unique tx values.                   #
        #--------------------------------------------------------------#
        tempTx <- data[use4fit,txName]
        levs <- sort(unique(round(tempTx,0L)))
        levs <- as.character(levs)
      }
      
      #------------------------------------------------------------------#
      # Reset treatment variable in data.frame using new vector.         #
      #------------------------------------------------------------------#
      tData[,txInfo@txName] <- tempTx
      
      #------------------------------------------------------------------#
      # obtain fit.                                                      #
      #------------------------------------------------------------------#
      fitResult <- Fit(object = moPropen[[k]], 
                       data = tData, 
                       response = tData[,txName])

      #------------------------------------------------------------------#
      # Store result as a new PropenSubsetFit object.                    #
      #------------------------------------------------------------------#
      propenFitObj[[k]] <- new("PropenSubsetFit",
                               subset = modelSubset,
                               levels = levs,
                               small = small,
                               modelObjectFit = fitResult)

    }

    #----------------------------------------------------------------------#
    # Recast fit objects to PropenSubsetFitList.                           #
    #----------------------------------------------------------------------#
    propenFitObj <- new("PropenSubsetFitList",
                        loo = propenFitObj,
                        txInfo = txInfo)

  } else if( is(moPropen, "modelObj") ){

    #----------------------------------------------------------------------#
    # Eliminate patients with only 1 tx option from dataset for fit        #
    #----------------------------------------------------------------------#
    use4fit <- eliminateSingleTx(subsets=subsets, ptsSubset=ptsSubset)

    #----------------------------------------------------------------------#
    # Determine if the smallest or largest treatment is not present        #
    # in returned predictions. Note that some methods return               #
    # predictions for all treatment options.                               #
    #----------------------------------------------------------------------#
    predArgs <- predictorArgs(moPropen)
    if( is(predArgs[["propen.missing"]], "NULL") ) {
      small <- TRUE
    } else if( tolower(predArgs[["propen.missing"]]) == "smallest" ) {
      small <- TRUE
      predArgs[["propen.missing"]] <- NULL
    } else if( tolower(predArgs[["propen.missing"]]) == "largest" ) {
      small <- FALSE
      predArgs[["propen.missing"]] <- NULL
    } 

    predictorArgs(moPropen) <- predArgs

    #----------------------------------------------------------------------#
    # feasibility rules also causes problems w/ identifiability of levels  #
    # The following tries to track what is being estimated. If tx is       #
    # a factor, the tx vector for the subset of data is recast as          #
    # another factor variable to get the correct levels. If tx is an       #
    # integer, the unique tx values for the subset of data are sorted      #
    # It is thus assumed that the prediction method will return the        #
    # probability matrix in the default order of factors or in sorted      #
    # order of the integer values.                                         #
    #----------------------------------------------------------------------#
    if( is.factor(data[,txName]) ) {
      #------------------------------------------------------------------#
      # If treatment is a factor, pull cases that are to be fit and      #
      # reset levels                                                     #
      #------------------------------------------------------------------#
      levs <- factor(data[use4fit,txName])
      levs <- levels(levs)
    } else {
      #------------------------------------------------------------------#
      # If treatment is not a factor, pull cases that are to be fit      #
      # and define levels as the unique values.                          #
      #------------------------------------------------------------------#
      levs <- sort(unique(round(data[use4fit,txName],0L)))
      levs <- as.character(levs)
    }


    #----------------------------------------------------------------------#
    # fit propen model using provided solver                               #
    #----------------------------------------------------------------------#
    propenFitObj <- fit(object = moPropen, 
                        data = data[use4fit,,drop=FALSE], 
                        response = data[use4fit,txName])

    #----------------------------------------------------------------------#
    # Store as PropenModelObjFit object.                                   #
    #----------------------------------------------------------------------#
    propenFitObj <- new("PropenModelObjFit",
                        modelObjectFit = propenFitObj,
                        levels = levs,
                        small = small)
  } else {

    DeveloperError("unrecognized class for moPropen", "fitPropen")

  }

  #--------------------------------------------------------------------------#
  # Recast result as one of class PropenFit.                                 #
  #--------------------------------------------------------------------------#
  result <- new("PropenFit",
                fits = propenFitObj)

  return(result)
}

