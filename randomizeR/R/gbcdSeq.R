#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class gbcdSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for gbcdSeq
# --------------------------------------------

# Representation of sequences for Generalized Biased Coin Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Generalized Biased Coin Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot rho a positive parameter which my be adjusted according to how strongly it is desired to balance the experiment.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
setClass("gbcdSeq", slots=c(rho = "numeric"), contains = "randSeq")


# --------------------------------------------
# Class definition for rgbcdSeq
# --------------------------------------------

# Representation of sequences for Efron's Biased Coin Design
# 
# @description This set of classes provides functionality of storing random 
# randomization sequences of Generalized Biased Coin Design along with the parameters 
# representing the design.
setClass("rGbcdSeq", contains = c("rRandSeq", "gbcdSeq"))

# --------------------------------------------
# Methods for gbcdSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "gbcdSeq"),
          function(obj) {
            if(obj@K == 2) {
              apply(obj@M, 1, function(x, rho){
                N <- length(x); p <- numeric(N); reps <- 1; p[1] <- 1/2
                while(reps < N){
                  nB <- sum(x[1:reps])
                  nA <- reps - nB
                  reps <- reps + 1
                  if(abs(nA-nB) == 0){
                    p[reps] <- 1/2
                  } else {
                    phi <- nB^rho/(nA^rho+nB^rho)
                    p[reps] <- (1-x[reps])*phi + x[reps]*(1-phi)
                  }
                }
                prod(p)
              }, rho <- obj@rho)
            } else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "gbcdSeq"),
          function(obj) {
            paste("gbcd(", round(obj@rho, digits = 2), ")", sep = "")
          }
)

