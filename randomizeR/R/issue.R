#' @include selBias.R
#' @include chronBias.R
#' @include bias.R
#' @include corGuess.R
#' @include imbalance.R
#' @include power.R
NULL

###############################################
# --------------------------------------------#
# Class issue                                 #
# --------------------------------------------#
###############################################

#' Issues in clinical trials
#' 
#' Summarizes the criteria for the assessment of randomization procedures.
#' 
#' @details
#' Randomization in clinical trials is supposed to 
#' control certain issues in clinical trials. 
#' Many of the issues are working in opposite direction, so it is crucial to decide for
#' which of the issues is relevant in the present clinical trial.
#' These issues
#' include
#' \itemize{
#' \item \strong{Selection bias} 
#' 		can occur if future treatment allocations are predictable due to 
#' 		restricted randomization and unmasking of past treatment assigments.
#'		Selection bias is represented by the \code{\link{selBias}} class.
#' \item \strong{Chronological bias} 
#' 		can occur if a time trend is present in the data. Time trends occur
#' 		due to learning curves, relaxed inclusion/ exclusion criteria or
#' new co-medication.
#'		Chronological bias is represented by the \code{\link{chronBias}} class.
#' \item \strong{Balance}
#' 		is important in order to ensure proper power estimation properties of
#' the treatments.
#' 		However, a high degree of balance favours selection bias.
#' A middle course seems optimal.
#'		Imbalance bias is represented by the \code{\link{imbal}} class. 
#' }
#' @name issue
#'
#' @seealso Representation of randomization procedures: \code{\link{randPar}}
#' @seealso Generation of randomization sequences: \code{\link{genSeq}}
#' @seealso Assessment of randomization sequences: \code{\link{assess}}
#' @seealso Comparison of randomization sequences: \code{\link{compare}}
#'
#' @family issues
#' 
#' @aliases issues
NULL

# Issue class
# 
# @name issue
setClassUnion("issue", c("selBias", "chronBias", "corGuess", "imbal", "power"))
#setClassUnion("issue", c("bias", "corGuess", "imbal", "power"))


# --------------------------------------------
# Accesssor functions for issue
# --------------------------------------------

#' Method defining the $ operator for the issue class
#' 
#' @inheritParams overview
setMethod("$", "issue",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for issue
# --------------------------------------------

setMethod("show", signature = "issue", definition = function(object) {
  validObject(object)
  # headline
  cat("\n Object of class \"", class(object)[1], "\"\n\n", sep = "")
  # iterate through all slots of the object
  names <- slotNames(object) 
  for(name in names){
    cat("\t", toupper(name), "=", slot(object, name), "\n")
  }
  cat("\n") 
})








