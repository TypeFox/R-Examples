#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# optimalClass_IPWE : calculates the IPWE contrast function for a single       #
#                      decision point binary tx.                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# txInfo   : an object of class txInfo                                         #
#                                                                              #
# propensity: A matrix of propensity scores.                                   #
#                                                                              #
# data     : data frame of covariates                                          #
#                                                                              #
# response : a response vector                                                 #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list                                                             =#
#=    constrast, mean.mu0                                                     =#
#=                                                                            =#
#==============================================================================#
optimalClass_IPWE <- function(txInfo, 
                              propensity, 
                              data,
                              response){

  #--------------------------------------------------------------------------#
  # Extract observed treatment                                               #
  #--------------------------------------------------------------------------#
  tx <- data[,TxName(txInfo)]

  #--------------------------------------------------------------------------#
  # Calculate IPWE contrast function.                                        #
  #--------------------------------------------------------------------------#
  ym <- tx/propensity[,"1"]*response -
        (1.0 - tx)/propensity[,"0"]*response

  #--------------------------------------------------------------------------#
  # Calculate non-contrast contribution to IPWE estimator.                   #
  #--------------------------------------------------------------------------#
  mmu <- (1.0 - tx)/propensity[,"0"]*response
  mmu <- sum(mmu)/nrow(data)

  return(list("contrast" = ym,
              "mean.mu0" = mmu))
}
