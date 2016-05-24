#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class crPar                                 #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the crPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatecrPar <- function(object) {
  errors <- character()
  
  if(length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for crPar
# --------------------------------------------

# Randomization paramters generic
setClass("crPar", contains = "randPar", validity = validatecrPar)


# --------------------------------------------
# Constructor function for crPar
# --------------------------------------------

#' Representing Complete Randomization
#' 
#' Represents the randomization procedure Complete Randomization.
#'
#' @details
#' Toss a fair coin \code{N} times in case \code{K=2} and assign the treatments according to the result of the coin. In case of \code{K>2}, replace the coin by a die with \code{K} sides.
#'
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @return 
#' \code{S4} object of the class \code{crPar}.
#'
#' @export
#'
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
crPar <- function(N, K = 2, ratio = rep(1, K), groups = LETTERS[1:K]) {
  new("crPar", N = N, K = K, ratio = ratio, groups = groups)
}


# --------------------------------------------
# Sampling algorithm for CR
# --------------------------------------------

# Complete Randomization
#
# This function implements a generalized version of Complete Randomization.
# In the original version Complete Randomization is equivalent to tossing
# a fair coin for a 1:1 allocation of subjects into two treatment groups. 
# This version extends the original version to support more than two treatment 
# groups and unequal allocation ratios.
#
# @inheritParams overview
# 
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain the number {j-1} at position \code{i}, 
# when patient \code{i} is allocated to treatment \code{j} (\code{j=1,...,K}).
completeRand <- function(N, K = 2, ratio = rep(1, K)) {
  sample(0:(K-1), size = N, prob = ratio/sum(ratio), replace = TRUE)
}


# --------------------------------------------
# Methods for crPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", 
          signature(obj = "crPar"),
          function(obj) {
            stopifnot(validObject(obj))
            if(obj@K != 2 || !identical(obj@ratio, c(1,1))) {
              stop("Only possible for K equals 2 and ratio corresponds to c(1,1).")
            }  
            new("crSeq", 
                M = compltSet(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "crPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
	    set.seed(seed)
            new("rCrSeq", 
                M = t(sapply(1:r, function(x) {
                  completeRand(N = N(obj), K = K(obj), ratio = ratio(obj))
                  })),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "crPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
			      set.seed(seed)
            new("rCrSeq", 
                M = t(completeRand(N = N(obj), K = K(obj), ratio = ratio(obj))), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "crPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rCrSeq", 
                M = t(sapply(1:r, function(x) {
                  completeRand(N = N(obj), K = K(obj), ratio = ratio(obj))
                  })),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "crPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rCrSeq", 
                M = t(completeRand(N = N(obj), K = K(obj), ratio = ratio(obj))), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "crPar"),
          function(obj) {
            "CR"
          }
)

