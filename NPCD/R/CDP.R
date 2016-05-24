
###############################################################################
# CDP:                                                                        #
#                                                                             #
# Compute probability of correct response of one item for one person given    #
# the item parameters, Q vector, and alpha vector.                            #
#                                                                             #
# Input:                                                                      #
# (1) Q: the Q-vector of the item. Columns represent attributes.              #
# (2) par: a list of parameters.                                              #
#          DINA&DINO --- par$slip: a scaler slip parameter for the item       #
#                   par$guess: a scaler guessing parameter for the item       #
#          NIDA --- par$slip: a vector of slip parameters for each attributes #
#                   par$guess: a vector of guess parameters for each attribute#
#          GNIDA --- par$slip: a vector of slip parameters for each attribute #
#                             for the item                                    #
#                    par$guess:  vector of guess parameters for each attribute#
#                             for the item                                    #
#          RRUM --- par$pi: a scaler pi parameter for the item                #
#                   par$r: a vector of r parameters for each attribute for the#
#                          item                                               #
# (3) alpha: a vector of attribute profile of the person.                     #
# (4) model: "DINA", "DINO", "NIDA", "GNIDA", "RRUM"                          #
#                                                                             #
# Output:                                                                     #
# (1) P: the probability of correct response for the item by the person       #
#                                                                             #
###############################################################################

CDP <- function(Q, par, alpha, model=c("DINA", "DINO", "NIDA", "GNIDA", "RRUM")){
  
  natt <- length(Q)
  model <- match.arg(model)
  
  if (model %in% c("DINA", "DINO", "NIDA", "GNIDA")){
	  par$slip[par$slip == 0] <- 0.001
	  par$guess[par$guess == 0] <- 0.001
	  par$slip[par$slip == 1] <- 0.999
	  par$guess[par$guess == 1] <- 0.999
  }
  
  if (model == "RRUM"){
	  par$pi[par$pi == 0] <- 0.001
	  par$r[par$r == 0] <- 0.001
	  par$pi[par$pi == 1] <- 0.999
	  par$r[par$r == 1] <- 0.999
  }
  
  {
  if (model == "DINA") 
  {
    ita <- prod(alpha ^ Q)
    P <- (1 - par$slip) ^ ita * par$guess ^ (1 - ita)
    
  } else if (model == "DINO")
  {
    omega <- 1 - prod((1 - alpha) ^ Q)
    P <- (1 - par$slip) ^ omega * par$guess ^ (1 - omega)  

  } else if (model %in% c("NIDA", "GNIDA"))
  {
    P <- prod(((1 - par$slip) ^ alpha * par$guess ^ (1 - alpha)) ^ Q)
      
  } else if (model == "RRUM")
  {
    P <- par$pi * prod(par$r ^ (Q * (1 - alpha)))
    
  } else return(warning("Model specification is not valid."))
 }
  
  return(P)
}
