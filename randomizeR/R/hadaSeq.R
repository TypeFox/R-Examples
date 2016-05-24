#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class hadaSeq                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for hadaSeq
# --------------------------------------------

# Representation of sequences for Hadamard Randomization
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Hadamard Randomization along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
setClass("hadaSeq", contains = "randSeq")


# --------------------------------------------
# Class definition for rHadaSeq
# --------------------------------------------

# Representation of sequences for Hadamard Randomization
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of Hadamard Randomization along with the parameters 
# representing the design.
setClass("rHadaSeq", contains = c("rRandSeq","hadaSeq"))


# --------------------------------------------
# Methods for hadaSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "hadaSeq"),
          function(obj){
            r <- ncol(obj@M)
            # every sequence equally probable, generate same probability
            rep(1/r, r)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "hadaSeq"),
          function(obj) {
            "HADA"
          }
)
