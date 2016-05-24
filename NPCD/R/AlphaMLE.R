
###############################################################################
# AlphaMLE                                                                    #
#                                                                             #
# Maximum likelihood estimation of attribute profiles for cogitive diagnostic #
# models.                                                                     #
#                                                                             #
# Input:                                                                      #
# (1) Y: a matrix of binary responses (1=correct, 0=incorrect). Rows          #
#               represent persons and columns represent items.                #
# (2) Q: the Q-matrix of the test. Rows represent items and columns represent #
#        attributes.                                                          #
# (3) par: a list of parameters.                                              #
#          DINA & DINO --- par$slip: a vector of slip parameters for each item#
#                   par$guess: a vector of guessing parameters for each item  #
#          NIDA --- par$slip: a vector of slip parameters for each attributes #
#                   par$guess: a vector of slip parameters for each attributes#
#          GNIDA --- par$slip: a matrix (item by attribute) of slip parameters#
#                    par$guess: a matrix (item by attribute) of guess param.  #
#          RRUM --- par$pi: a vector of the pi parameters for each item       #
#                   par$r: a matrix (# items by # attributes) of r parameters #
# (4) model: "DINA", "DINO", "NIDA", "GNIDA", "RRUM"                          #
# (5) NP.method: "Hamming", "Weighted", "Penalized"                           #
# (6) undefined.flag: a binary vector indicating whether the parameters of    #
#                     each item are undefined (1=undefined, 0=defined).       #
#                                                                             #
# Output:                                                                     #
# (1) alpha.est: a matrix of estimated attribute profiles for all examinees   #
# (2) pattern: all attribute profiles in the search space.                    #
# (3) loglike.matrix: The values for the loss function. Rows represent        #
#                     candidate attribute profiles in the same order with the #
#                     pattern matrix; Columns represent different examinees.  #
#                                                                             #
###############################################################################

AlphaMLE <- function(Y, Q, par, model=c("DINA", "DINO", "NIDA", "GNIDA", "RRUM"), undefined.flag=NULL) {
  
  #####
  # 1 #
  ##### Check dimension consistency and convert data to the right formats 
  
  Y <- as.matrix(Y)
  Q <- as.matrix(Q) 
  check <- NULL
  check <- CheckInput(Y, Q)  
  if (!is.null(check)) return(warning(check))
  
  model <- match.arg(model)
  
  #####
  # 2 #
  ##### Estimation
  
  nperson <- dim(Y)[1]
  nitem <- dim(Q)[1]
  natt <- dim(Q)[2]
  pattern <-AlphaPermute(natt)
  loglike.matrix <- matrix(NA, dim(pattern)[1], nperson)  
  alpha.est <- matrix(NA, nperson, natt)
  est.class <- rep(NA, nperson)
  n.tie <- rep(0, nperson)
  class.tie <- matrix(0, nperson, dim(pattern)[1])
  
  for (i in 1:nperson)
  {
    loglike <- NULL
    
    for (j in 1:nrow(pattern))
    {
      loglike[j] <- CDL(Y[i, ], Q, par, pattern[j, ], model, undefined.flag)  
    }
    
    loglike.matrix[, i] <- loglike
    
    if (length(which(loglike == max(loglike))) == 1){
      est.class[i] <- which(loglike == max(loglike))
    } else {
      n.tie[i] <- length(which(loglike == max(loglike)))
      class.tie[i, 1:sum(loglike == max(loglike))] <- which(loglike == max(loglike))
      est.class[i] <- sample(which(loglike == max(loglike)), 1)      
    }
    alpha.est[i, ] <- pattern[est.class[i], ]
  }  
  
  output <- list(alpha.est=alpha.est, est.class=est.class, pattern=pattern, n.tie=n.tie, class.tie=class.tie, loglike.matrix=loglike.matrix, Y=Y, Q=Q, par=par, model=model)
  class(output) <- "AlphaMLE"
  return(output)
}
