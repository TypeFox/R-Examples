#' @include randPar.R
#' @include pbrPar.R
NULL

###############################################
# --------------------------------------------#
# Class rarPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the rarPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validaterarPar <- function(object) {
  errors <- character()
  N <- object@N
  ratio <- object@ratio
  
  if (!(N %% sum(ratio) == 0)) {
    msg <- paste("N = ", N, " is not a multiple of sum(ratio) = "
                 , sum(ratio),".", sep = "")
    errors <- c(errors, msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for rarPar
# --------------------------------------------

# Randomization paramters generic
setClass("rarPar", contains = "randPar", validity = validaterarPar)


# --------------------------------------------
# Constructor function for rarPar
# --------------------------------------------

#' Representing Random Allocation Rule
#' 
#' Represents the randomization procedure Random Allocation Rule.
#'
#' @details
#' Fix a total sample size \code{N} the number of treatment groups \code{K},
#' and the vector of the \code{ratio}. Afterwards, all patients are assigned
#' according to the ratio to the corresponding treatment groups.
#' All randomization sequences are equiprobable.
#'
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @return 
#' \code{S4} object of the class \code{rarPar}.
#'
#' @export
#'
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
rarPar <- function(N, K = 2, ratio = rep(1, K), groups = LETTERS[1:K]) {
  new("rarPar", N = N, K = K, ratio = ratio, groups = groups)
}


# --------------------------------------------
# Methods for rarPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", 
          signature(obj = "rarPar"),
          function(obj) {
            if(obj@K != 2 || !identical(obj@ratio, c(1,1))) {
              stop("Only possible for K equals 2 and ratio corresponds to c(1,1).")
            }  
            allSeqs <- compltSet(obj)
            finBal <- apply(allSeqs, 1, function(x) 2*sum(x) == length(x))
            new("rarSeq", 
                M = allSeqs[finBal, ] , 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "rarPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rRarSeq", 
                M = t(sapply(1:r, function(x) {
                  blockRand(bc = N(obj), K = K(obj), ratio = ratio(obj))
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
          signature(obj = "rarPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
	    set.seed(seed)
            new("rRarSeq", 
                M = t(blockRand(bc = N(obj), K = K(obj), ratio = ratio(obj))), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		        seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rarPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed = sample(1:(10^6),1)
            set.seed(seed)
            new("rRarSeq", 
                M = t(sapply(1:r, function(x) {
                  blockRand(bc = N(obj), K = K(obj), ratio = ratio(obj))
                  })), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		        seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rarPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed = sample(1:(10^6),1)
	    set.seed(seed)
            new("rRarSeq", 
                M = t(blockRand(bc = N(obj), K = K(obj), ratio = ratio(obj))), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		        seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "rarPar"),
          function(obj) {
            "RAR"
          }
)
