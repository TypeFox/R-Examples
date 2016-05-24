#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class bbcdSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for bbcdSeq
# --------------------------------------------

# Representation of sequences for Bayesian Biased Coin Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Bayesian Biased Coin Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot a a positive parameter which my be adjusted according to how strongly it is desired to balance the experiment.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
setClass("bbcdSeq", slots=c(a = "numeric"), contains = "randSeq")


# --------------------------------------------
# Class definition for rbbcdSeq
# --------------------------------------------

# Representation of sequences for Bayesian Biased Coin Design
# 
# @description This set of classes provides functionality of storing random 
# randomization sequences of Bayesian Biased Coin Design along with the parameters 
# representing the design.
setClass("rbbcdSeq", contains = c("rRandSeq", "bbcdSeq"))

# --------------------------------------------
# Methods for bbcdSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "bbcdSeq"),
          function(obj) {
            if(obj@K == 2) {
              N <- ncol(obj@M)
              apply(obj@M, 1, function(x, a){
                p <- numeric(N); p[1] <- 1/2; p[2] <- 1; reps <- 2
                while(reps < N){
                  nB <- sum(x[1:reps])
                  nA <- reps - nB
                  if (nA == nB) {
                    p[reps + 1] <- 1/2 
                  } else {
                    f <- (1+nB/(reps*nA))^(1/a)/((1+nB/(reps*nA))^(1/a) + (1+nA/(reps*nB))^(1/a))
                    p[reps + 1] <- (1-x[reps+1])*f + x[reps+1]*(1-f)
                  }
                  reps <- reps + 1 
                }
                prod(p)
              }, a <- obj@a)
            } else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "bbcdSeq"),
          function(obj) {
            paste("bbcd(", round(obj@a, digits = 2), ")", sep = "")
          }
)

