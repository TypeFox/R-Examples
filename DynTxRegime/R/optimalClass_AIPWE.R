#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# optimalClass_AIPWE : calculates the AIPWE contrast function for a single     #
#                      decision point binary tx.                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# outcome   : an object of type SimpleFit or IterateFit                        #
#                                                                              #
# txInfo    : an object of class TxInfo                                        #
#                                                                              #
# propensity: A matrix of propensity scores.                                   #
#                                                                              #
# data      : data frame of covariates                                         #
#                                                                              #
# response  : a response vector                                                #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list                                                             =#
#=    constrast, mean.mu0                                                     =#
#=                                                                            =#
#==============================================================================#
optimalClass_AIPWE <- function(outcome, 
                               txInfo, 
                               propensity, 
                               data,
                               response){

  #--------------------------------------------------------------------------#
  # Extract treatment options                                                #
  #--------------------------------------------------------------------------#
  sset <- SuperSet(txInfo)

  #--------------------------------------------------------------------------#
  # Extract observed treatment                                               #
  #--------------------------------------------------------------------------#
  tx <- data[,TxName(txInfo)]

  n <- nrow(data)

  #--------------------------------------------------------------------------#
  # Duplicate dataset for each treatment option  {0,1}                       #
  #--------------------------------------------------------------------------#
  dft <- rbind(data, data)
  dft[,TxName(txInfo)] <- c(rep(0L,n), rep(1L,n))

  #--------------------------------------------------------------------------#
  # Predict outcome                                                          #
  #--------------------------------------------------------------------------#
  me <- PredictMain(object=outcome, newdata=dft)
  cn <- PredictCont(object=outcome, newdata=dft)

  #--------------------------------------------------------------------------#
  # Recast as a two column matrix; one column for each treatment.            #
  #--------------------------------------------------------------------------#
  mu <- matrix(me + cn, ncol = 2L)

  #--------------------------------------------------------------------------#
  # Calculate AIPWE contrast function.                                       #
  #--------------------------------------------------------------------------#
  ym <- tx/propensity[,"1"]*response -
        (1.0 - tx)/propensity[,"0"]*response -
        (tx - propensity[,"1"])/propensity[,"1"]*mu[,2L] - 
        (tx - propensity[,"1"])/propensity[,"0"]*mu[,1L]

  #--------------------------------------------------------------------------#
  # Calculate non-contrast contribution to AIPWE estimator.                  #
  #--------------------------------------------------------------------------#
  mmu <- (1.0 - tx)/propensity[,"0"]*response +
         (tx - propensity[,"1"])/propensity[,"0"]*mu[,1L]

  mmu <- sum(mmu)/n

  return(list("contrast" = ym,
              "mean.mu0" = mmu))
}
