#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# iqLearn_pm : Calculate value functions for tx = +1/-1 for new patient        #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# object : object of class iqLearnSS, iqLearnFS_ME, iqLearnFS_C, or            #
#          iqLearnFS_VHet                                                      #
#                                                                              #
# newdata : data.frame of covariate information for new patient                #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= returns a matrix of value functions column 1: tx=-1, column 2: tx = 1.     =#
#=                                                                            =#
#==============================================================================#
iqLearn_pm <- function(object, 
                       newdata){

  #--------------------------------------------------------------------------#
  # Retrieve the number of samples in new dataset.                           #
  #--------------------------------------------------------------------------#
  n <- nrow(newdata)

  #--------------------------------------------------------------------------#
  # Replicate data to calculate +/- simultaneously.                          #
  #--------------------------------------------------------------------------#
  data <- rbind(newdata, newdata)

  #--------------------------------------------------------------------------#
  # First n samples have tx = -1. Second n samples have tx = +1              #
  #--------------------------------------------------------------------------#
  tx <- c(rep(-1L,n),rep(1L,n))

  if( TxName(object) %in% colnames(newdata) ){
    data[,TxName(object)] <- tx
  } else {
    data <- cbind(data,tx)
    colnames(data) <- c(colnames(newdata),TxName(object))
  }

  #--------------------------------------------------------------------------#
  # Calculate value function                                                 #
  # Note that predict method has to be called because combined fits cannot   #
  # be broken down into "main" and "contrast"; "contrast" will be sent back  #
  # as zero.                                                                 #
  # For iterate fits, the formula has been modified to include the treatment #
  # variable, therefore the value of treatment does not need to be considered#
  #--------------------------------------------------------------------------#
  mn <- PredictMain(object = object, newdata = data)
  cn <- PredictCont(object = object, newdata = data)

  #--------------------------------------------------------------------------#
  # Return positive and negative value functions                             #
  #--------------------------------------------------------------------------#
  return( matrix(data = mn + cn, ncol = 2L) )
}

