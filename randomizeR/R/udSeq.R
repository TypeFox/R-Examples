#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class udSeq                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for udSeq
# --------------------------------------------

# Representation of sequences for Wei's Urn Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Wei's Urn Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot ini integer representing the initial urn composition.
# @slot add integer representing the number of balls that are added to the
# urn in each step.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
#
setClass("udSeq", slots= c(ini = "numeric", add = "numeric"), contains = "randSeq")


# --------------------------------------------
# Class definition for rUdSeq
# --------------------------------------------

# Representation of sequences for Wei's Urn Design
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of Wei's Urn Design along with the parameters 
# representing the design.
setClass("rUdSeq", contains = c("rRandSeq", "udSeq"))

# --------------------------------------------
# Methods for udSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "udSeq"),
          function(obj) {
            if(obj@K == 2) {
              apply(obj@M, 1, function(randSeq, ini, add) {
                conditionalProb <- numeric(length(randSeq))
                conditionalProb[1] <- 1/2
                sumR <- randSeq[1]
                for(j in 1:(length(randSeq) - 1)) {
                  biasedCoin <- (ini + add*(j-sumR))/(2*ini + add*j)
                  if(randSeq[j+1] == 1) {
                    conditionalProb[j+1] <- biasedCoin
                  } else {
                    conditionalProb[j+1] <- 1 - biasedCoin
                  }
                  sumR <- sumR + randSeq[j+1]
                }
                return(prod(conditionalProb))
              }, ini = obj@ini, add = obj@add)
            } 
            else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "udSeq"),
          function(obj) {
            paste("UD(", obj@ini, ",", obj@add, ")", sep = "")
          }
)
