#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class abcdPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the abcdPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateabcdPar <- function(object) {
  errors <- character()
  a <- object@a
  lengtha <- length(a)
  
  if(lengtha != 1) {
    msg <- paste("a has length ", lengtha, ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if (a < 0) {
    msg <- paste("First element of a is ", a, ". Should be positive.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for abcdPar
# --------------------------------------------

# Randomization paramters generic
setClass("abcdPar",
         slots = c(a = "numeric"),
         contains = "randPar",
         validity = validateabcdPar)

# --------------------------------------------
# Constructor function for abcdPar
# --------------------------------------------

#' Representing Adjustable Biased Coin Design
#' 
#' Represents the randomization procedure Adjustable Biased Coin Design.
#'
#' @details
#' This is a class of 'biased coins' where the probability of selecting the under-represented
#' treatment is dependent from the absolute difference between the two treatment allocations
#' up to the current step. 
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @references
#' A. B. Antognini and A. Giovagnoli (2004) A new 'biased coin design' for the sequential 
#' allocation of two treatments.
#' \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}
#' \strong{53}, No. 4, 651-664
#' 
#' @return 
#' \code{S4} object of the class \code{abcdPar}.
#'
#' @export
abcdPar <- function(N, a, groups = LETTERS[1:2]) {
  new("abcdPar", N = N, a = a, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Sampling algorithm for abcd
# --------------------------------------------

#' Sampling algorithm for abcd 
#'
#' @inheritParams overview
#' 
#' @return A vector with the allocation sequence for a clinical trial. 
#' It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
#' is allocated to treatment A (resp. B).
#' 
#' @export
#' 
#' @references
#' A. B. Antognini and A. Giovagnoli (2004) A new 'biased coin design' for the sequential 
#' allocation of two treatments.
#' \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}
#' \strong{53}, No. 4, 651-664
#' 
abcdRand <- function(N, a, K = 2) {
  stopifnot(round(N) == N)
  if (!(K == 2)) stop("abcd: K > 2 not available yet.")
  R <- numeric(N); reps <- 0; nA <- 0; nB <- 0; d <- 0
  while(reps < N) {
    # case analysis
    if (abs(d) == 0) {
      R[reps + 1] <- rbinom(1, 1, 0.5)  
    } else {
      if(d>=1){
        Fa <- 1/(abs(d)^a + 1)
      } else {
        Fa <- (abs(d)^a)/(abs(d)^a + 1)
      }
      R[reps + 1] <- rbinom(1, 1, 1 - Fa)  
    }
    reps <- reps + 1
    nB <- sum(R[1:reps])
    nA <- reps - nB
    d <- nA-nB
  }
  R
}



# --------------------------------------------
# Methods for abcdPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "abcdPar"),
          function(obj) {
            new("abcdSeq", 
                M = compltSet(obj),
                a = a(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups
            )
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "abcdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rAbcdSeq", 
                M = t(sapply(1:r, function(x) {
                  abcdRand(N = N(obj), a = a(obj))
                })),
                N = N(obj),
                a = a(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
                seed = seed)
          }
)
#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "abcdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rAbcdSeq", 
                M = t(sapply(1:r, function(x) {
                  abcdRand(N = N(obj), a = a(obj))
                })),
                N = N(obj),
                a = a(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "abcdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rAbcdSeq", 
                M = t(abcdRand(N = N(obj), a = a(obj))), 
                a = a(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "abcdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rAbcdSeq", 
                M = t(abcdRand(N = N(obj), a = a(obj))), 
                a = a(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "abcdPar"),
          function(obj) {
            paste("abcd(", obj@a, ")", sep = "")
          }
)

