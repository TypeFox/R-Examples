#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class abcdSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for abcdSeq
# --------------------------------------------

# Representation of sequences for Adjustable Biased Coin Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Adjustable Biased Coin Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot a a positive parameter which my be adjusted according to how strongly it is desired to balance the experiment.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
setClass("abcdSeq", slots=c(a = "numeric"), contains = "randSeq")


# --------------------------------------------
# Class definition for rabcdSeq
# --------------------------------------------

# Representation of sequences for Adjustable Biased Coin Design
# 
# @description This set of classes provides functionality of storing random 
# randomization sequences of Adjustable Biased Coin Design along with the parameters 
# representing the design.
setClass("rAbcdSeq", contains = c("rRandSeq", "abcdSeq"))

# --------------------------------------------
# Methods for abcdSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "abcdSeq"),
          function(obj) {
            if(obj@K == 2) {
              apply(obj@M, 1, function(x, a){
                N <- length(x); p <- numeric(N); reps <- 1; p[1] <- 1/2
                while(reps < N){
                  nB <- sum(x[1:reps])
                  nA <- reps - nB
                  d <- nA-nB
                  reps <- reps + 1
                  if(abs(d) == 0){
                    p[reps] <- 1/2
                  } else {
                    if(d>=1){
                      Fa <- 1/(abs(d)^a + 1)
                    } else {
                      Fa <- (abs(d)^a)/(abs(d)^a + 1)
                    }
                    p[reps] <- (1-x[reps])*Fa + x[reps]*(1 - Fa)
                  }
                }
                prod(p)
              }, a <- obj@a)
            } else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "abcdSeq"),
          function(obj) {
            paste("abcd(", round(obj@a, digits = 2), ")", sep = "")
          }
)

