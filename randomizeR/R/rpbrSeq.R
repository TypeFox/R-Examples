#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class rpbrSeq                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for rpbrSeq
# --------------------------------------------

# Representation of sequences for Permuted Block Randomization
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Permuted Block Randomization along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot bc vector which contains the lengths \code{k_1,...,k_l} of each block. 
# @slot M matrix containing randomisation sequences of length \code{N} in its rows.
#
setClass("rRpbrSeq", slots = c(rb = "numeric", bc = "list", filledBlock = "logical"),  
         contains = "rRandSeq")


# --------------------------------------------
# Methods for rpbrSeq
# --------------------------------------------

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "rRpbrSeq"),
          function(obj) {
            rb <- capture.output(cat(obj@rb, sep = ","))
            if (obj@filledBlock) {
              paste(c("RPBRFB(", rb, ")"), sep = "", collapse = "")
            } else {
              paste(c("RPBR(", rb, ")"), sep = "", collapse = "")
            }  
          }
)
