#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class ebcPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the ebcPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateebcPar <- function(object) {
  errors <- character()
  p <- object@p
  N <- object@N
  ratio <- object@ratio

  if(p[1] < 0.5 || p[1] > 1) {
    msg <- paste("First element of p is ", p[1], ". Should be in [0.5,1].", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(length(p) > 1) {
    msg <- paste("p has length  ", length(p), ". Should be one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for ebcPar
# --------------------------------------------

# Randomization paramters generic
setClass("ebcPar",
         slots = c(p = "numeric"),
         contains = "randPar",
         validity = validateebcPar)

# --------------------------------------------
# Constructor function for ebcPar
# --------------------------------------------

#' Representing Efron's Biased Coin Design
#' 
#' Represents the randomization procedure Efron's Biased Coin Design.
#'
#' @details
#' Flip a biased coin with probability \code{p} in favour of the treatment
#' which is allocated less frequently. If both treatments have been assigned
#' equally often a fair coin is tossed.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#' 
#' @return 
#' \code{S4} object of the class \code{ebcPar}.
#'
#' @references 
#' B. Efron (1971) Forcing a sequential experiment to be balanced. \emph{Biometrika},
#' \strong{58}, 403-17.
#' 
#' @export
ebcPar <- function(N, p, groups = LETTERS[1:2]) {
  new("ebcPar", N = N, p = p, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Sampling algorithm for EBC
# --------------------------------------------

# Efrons Biased Coin and Big Stick Design 
#
# This procedure generalises efrons biased coin design. It permits a maximum
# tolerated imbalance \code{MTI} during the trial.
# In the setting with success probability p = 0.5 of the biased coin it thus
# yields the Big Stick Design.
#
# @inheritParams overview
# 
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
# is allocated to treatment A (resp. B).
# 
# @export
# 
# @references 
# B. Efron (1971) Forcing a sequential experiment to be balanced. \emph{Biometrika},
# \strong{58}, 403-17.
# J. F. Soares and C. F. Jeff Wu (1983) Some Restricted Randomization Rules in
# Sequential Designs. \emph{Comm. in Stat.}, \strong{12}, 2017-34. 
efronRand <- function(bc, p, mti, K = 2) {
  stopifnot(is.numeric(bc), is.numeric(p), round(mti) == mti, p >= 0.5 & p <= 1)
  if(!(K == 2)) stop("EBC, BSD: K>2 not available yet.")
  N <- sum(bc)
  R <- numeric(N); reps <- 0; sumR <- 0; imb <- 0
  while(reps < N) {
    # case analysis
    if(imb == 0) R[reps + 1] <- rbinom(1, 1, 0.5)  
    else if (imb >= mti) R[reps + 1] <- 0
    else if (imb <= -mti) R[reps + 1] <- 1  
    else if (sign(imb) == 1) R[reps + 1] <-  rbinom(1, 1, 1 - p)
    else R[reps+1] <- rbinom(1, 1, p)
    
    reps <- reps + 1
    sumR <- sumR + R[reps]
    imb <- 2*sumR - reps
  }
  R
}


# --------------------------------------------
# Methods for ebcPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "ebcPar"),
          function(obj) {
            new("ebcSeq", M = compltSet(obj), p = coin(obj),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "ebcPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rEbcSeq", 
                M = t(sapply(1:r, function(x) {
                  efronRand(bc = N(obj), p = coin(obj), mti = N(obj), K = K(obj))
                          })),
                N = N(obj),
                p = coin(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
	              seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "ebcPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rEbcSeq", 
                M = t(efronRand(bc = N(obj), p = coin(obj), mti = N(obj))), 
                p = coin(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "ebcPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rEbcSeq", 
                M = t(sapply(1:r, function(x) {
                  efronRand(bc = N(obj), p = coin(obj), mti = N(obj), K = K(obj))
                          })),
                N = N(obj),
                p = coin(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "ebcPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rEbcSeq", 
                M = t(efronRand(bc = N(obj), p = coin(obj), mti = N(obj))), 
                p = coin(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "ebcPar"),
          function(obj) {
            paste("EBC(", round(obj@p, digits = 2), ")", sep = "")
          }
)



