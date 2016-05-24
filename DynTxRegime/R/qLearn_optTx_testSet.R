#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# qLearn_optTx_testSet :  Function to estimate optimal stage-n treatment rule  #
# for new patients based on information gleaned from training set              #
# This is a private function used by the public S4 method 'optTx'              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# object  : object of class QLearn created by a call to function qLearn        #
#                                                                              #
# newdata : covariates of new patient                                          #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
qLearn_optTx_testSet <- function (object, 
                                  newdata){

  n <- nrow(newdata)
  p <- ncol(newdata)

  #--------------------------------------------------------------------------#
  # Retrieve superset of tx options and tx column header from QLearn object  #
  #--------------------------------------------------------------------------#
  sset <- SuperSet(object)
  txN <- TxName(object)

  #--------------------------------------------------------------------------#
  # If tx column header not in newdata, add the column to data               #
  #--------------------------------------------------------------------------#
  if( !(txN %in% colnames(newdata)) ) {
    nms <- colnames(newdata)
    newdata <- cbind(newdata,c(rep(1L,n)))
    colnames(newdata) <- c(nms, txN)
  }

  #--------------------------------------------------------------------------#
  # Calculate Q-Functions are all tx options                                 #
  #--------------------------------------------------------------------------#
  qFunctions <- matrix(data = 0.0,
                       nrow = n,
                       ncol = length(sset),
                       dimnames = list(NULL, sset))

  for( i in 1:length(sset) ){
    if( is(Base(object),"character") ) {
      newdata[,txN] <- factor(sset[i], levels = sset)
    } else {
      newdata[,txN] <- as.integer(sset[i])
    }

    me <- PredictMain(object = object,
                      newdata = newdata)

    cn <- PredictCont(object = object,
                      newdata = newdata)

    qFunctions[,i] <- me + cn

  }

  #--------------------------------------------------------------------------#
  # Execute feasibility rules to determine what txs are available to each pt #
  #--------------------------------------------------------------------------#
  txList<- feasibility(superSet = SuperSet(object), 
                       fSet = SubsetRule(object), 
                       txName = TxName(object),
                       data = newdata)

  for( i in 1:length(txList$subsets) ) {

    #----------------------------------------------------------------------#
    # Determine which patients fall into this tx subset                    #
    #----------------------------------------------------------------------#
    tst <- txList$ptsSubset == names(txList$subsets)[i]

    #----------------------------------------------------------------------#
    # Extract tx values for this tx subset                                 #
    #----------------------------------------------------------------------#
    tx <- txList$subsets[[i]]

    #----------------------------------------------------------------------#
    # Identify the elements of the superset contained in the subset        #
    #----------------------------------------------------------------------#
    j <- sset %in% tx

    #----------------------------------------------------------------------#
    # Reset Q-Functions for those not in the subset to NA                  #
    #----------------------------------------------------------------------#
    qFunctions[tst,!j] <- NA

  }

  #--------------------------------------------------------------------------#
  # Identify the maximum Q-Function, thus optimal tx                         #
  #--------------------------------------------------------------------------#
  fun <- function(x){
    if(all(is.na(x))) { 
      return(NA) 
    } else {
      return(which.max(x))
    }
  }

  q2opt <- apply(qFunctions,1,fun)

  if( is(Base(object), "character") ) {
    optTx <- factor(colnames(qFunctions)[q2opt],
                    levels = colnames(qFunctions))
  } else {
    optTx <- as.integer(colnames(qFunctions)[q2opt])
  }

  #--------------------------------------------------------------------------#
  # Return all Q-functions and the treatment that yields the largest         #
  #--------------------------------------------------------------------------#
  return(list("qFunctions" = qFunctions,
               "optimalTx" = optTx))
}
