#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class pbrSeq                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for pbrSeq
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
setClass("pbrSeq", slots = c(bc = "numeric"),  contains = "randSeq")

# --------------------------------------------
# Class definition for rPbrSeq
# --------------------------------------------

# Representation of sequences for Permuted Block Randomization
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of Permuted Block Randomization along with the parameters 
# representing the design.
setClass("rPbrSeq",  contains = c("rRandSeq","pbrSeq"))


# --------------------------------------------
# Methods for pbrSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "pbrSeq"),
          function(obj) {
            if(obj@K == 2) {
             # every sequence equally probable, generate same probability
             r <- nrow(obj@M)
             bc <- blocks(obj)
             rep(1/prod(sapply(bc, function(x) choose(x, x/2))), r)
            }
            else "Only supported for K=2."
          } 
)


#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "pbrSeq"),
          function(obj) {
            if (sum(!duplicated(obj@bc)) == 1) {
              paste("PBR(", obj@bc[1], ")", sep = "")
            } else {
              bc <- capture.output(cat(obj@bc, sep = ","))
              paste(c("PBR(", bc, ")"), sep = "", collapse = "")
            }  
          }
)
