#' @include randSeq.R
#' @include util.R
NULL

###############################################
# --------------------------------------------#
# Class tbdSeq                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for tbdSeq
# --------------------------------------------

# Representation of sequences for the Truncated Binomial Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of the Truncated Binomial Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot bc vector which contains the realized block lengths from the trial.
# @slot M matrix containing randomisation sequences of length \code{N} in its
# rows.
setClass("tbdSeq", slots=c(bc = "numeric"), contains = "randSeq")

# --------------------------------------------
# Class definition for rTbdSeq
# --------------------------------------------
#
# Representation of sequences for the Truncated Binomial Design
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of the Truncated Binomial Design along with the parameters 
# representing the design.
setClass("rTbdSeq", contains = c("rRandSeq", "tbdSeq"))

# --------------------------------------------
# Methods for tbdSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "tbdSeq"),
          function(obj) {
            # every sequence equally probable, generate same probability
            seqs <- obj@M
            bc <- blocks(obj)
            
            apply(seqs, 1, function(x, bc){
              chunks <- split(x, rep(1:length(bc), bc))
              probs <- sapply(chunks, function(y) 0.5^lastRandom(y))
              prod(probs)
            }, bc)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "tbdSeq"),
          function(obj) {
            if (sum(!duplicated(obj@bc)) == 1) {
              paste("TBD(", obj@bc[1], ")", sep = "")
            } else {
              bc <- capture.output(cat(obj@bc, sep = ","))
              paste(c("TBD(", bc, ")"), sep = "", collapse = "")
            }  
          }
)


