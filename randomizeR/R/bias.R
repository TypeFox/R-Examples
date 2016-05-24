#' @include selBias.R
#' @include chronBias.R
NULL

###############################################
# --------------------------------------------#
# Class bias                                  #
# --------------------------------------------#
###############################################

# Bias class
# 
# @name bias
setClassUnion("bias", c("chronBias", "selBias"))


# --------------------------------------------
# Accesssor functions for bias
# --------------------------------------------

#' Get type of an object
#'
#' Accesses the type slot of an S4 object
#' 
#' @param obj a bias object (i.e. S4 object inheriting from \code{bias})
#' 
#' @return 
#' Character string specifying the type of bias \code{obj} represents, e.g. \code{"linT"} in 
#' case of chronological bias.
#'
#' @export
type <- function(obj) {
  if (.hasSlot(obj, "type")) obj@type else stop("Object has no slot named type.")   
} 



