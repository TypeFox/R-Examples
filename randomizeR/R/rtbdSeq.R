#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class rRtbdSeq                              #
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
# @slot bc vector which contains the realized block lengths from the trial.
# @slot p success probability of the coin.
# @slot rb random blocks to choose from
# 
setClass("rRtbdSeq", slots=c(rb = "numeric", bc = "list", filledBlock = "logical"),
         contains = "rRandSeq")

# --------------------------------------------
# Methods for rtbdSeq
# --------------------------------------------

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "rRtbdSeq"),
          function(obj) {
            rb <- capture.output(cat(obj@rb, sep = ","))
            if (obj@filledBlock) {
              paste(c("RTBDFB(", rb, ")"), sep = "", collapse = "")
            } else {
              paste(c("RTBD(", rb, ")"), sep = "", collapse = "")
            }
          }
)
