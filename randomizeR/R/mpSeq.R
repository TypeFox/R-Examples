#' @include randSeq.R
#' @include mpPar.R
NULL

###############################################
# --------------------------------------------#
# Class mpSeq                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for mpSeq
# --------------------------------------------

# Representation of sequences for the Maximal Procedure
# 
# @description This set of classes provides functionality of storing randomization
# sequences of the Maximal Procedure along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial
# @slot mti The maximum tolerated imbalance during the trial
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
# 
setClass("mpSeq", slots = c(mti = "numeric"), contains = "randSeq")

# --------------------------------------------
# Class definition for rMpSeq
# --------------------------------------------

# Representation of sequences for the Maximal Procedure
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of the Maximal Procedure along with the parameters 
# representing the design.
setClass("rMpSeq", contains = c("rRandSeq","mpSeq"))

# --------------------------------------------
# Methods for mpSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "mpSeq"),
          function(obj) {
            # every sequence euqiprobable, generate same probability
            r <- nrow(obj@M)
            tot <- (createMPMatrix(N(obj), mti(obj), ratio(obj)))[1]
            rep(1/tot, r)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "mpSeq"),
          function(obj) {
            paste("MP(", obj@mti, ")", sep = "")
          }
)

