#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class bbcdPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the bbcdPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatebbcdPar <- function(object) {
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
# Class definition for bbcdPar
# --------------------------------------------

# Randomization paramters generic
setClass("bbcdPar",
         slots = c(a = "numeric"),
         contains = "randPar",
         validity = validatebbcdPar)

# --------------------------------------------
# Constructor function for bbcdPar
# --------------------------------------------

#' Representing Bayesian Biased Coin Design
#' 
#' Represents the randomization procedure Bayesian Biased Coin Design.
#'
#' @details
#' Extension of Efron's biased coin design.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @references
#' A. B. Antognini and Maroussa Zagoraiou (2014) Balance and randomness in sequential
#' clinical trials: the dominant biased coin design.
#' \emph{Pharmaceutical Statistics}
#' \strong{13(2)}, 119-127
#' 
#' @return 
#' \code{S4} object of the class \code{bbcdPar}.
#'
#' @export
bbcdPar <- function(N, a, groups = LETTERS[1:2]) {
  new("bbcdPar", N = N, a = a, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Sampling algorithm for bbcd
# --------------------------------------------

#' Sampling algorithm for bbcd 
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
#' A. B. Antognini and Maroussa Zagoraiou (2014) Balance and randomness in sequential
#' clinical trials: the dominant biased coin design.
#' \emph{Pharmaceutical Statistics}
#' \strong{13(2)}, 119-127
#' 
bbcdRand <- function(N, a, K = 2) {
  stopifnot(round(N) == N)
  if (!(K == 2)) stop("bbcd: K > 2 not available yet.")
  R <- numeric(N); nA <- 0; nB <- 0
  
  # special case for n = 1 and n = 2 (PBR for n = 2)
  R[1] <- rbinom(1, 1, 0.5)
  nB <- sum(R[1])
  nA <- 1 - nB
  R[2] <- ifelse(nB == 1, 0, 1)
  reps <- 2
  
  # n > 2
  while(reps < N) {
    nB <- sum(R[1:reps])
    nA <- reps - nB
    # case analysis
    if (nA == nB) {
      R[reps+1] <- rbinom(1, 1, 0.5)  
    } else {
      f <- (1+nB/(reps*nA))^(1/a)/((1+nB/(reps*nA))^(1/a) + (1+nA/(reps*nB))^(1/a))
      R[reps + 1] <- rbinom(1, 1, 1 - f)  
    }
    reps <- reps + 1
  }
  R
}



# --------------------------------------------
# Methods for bbcdPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "bbcdPar"),
          function(obj) {
            allSeqs <- compltSet(obj) 
            # Two successive treatment assignments at the beginning of the trial 
            # cannot be the same 
            inside <- which(apply(allSeqs, 1, function(x) sum(x[1:2]) == 1) == TRUE)
            new("bbcdSeq", 
                M = allSeqs[inside,],
                a = a(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups
            )
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "bbcdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rbbcdSeq", 
                M = t(sapply(1:r, function(x) {
                  bbcdRand(N = N(obj), a = a(obj))
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
setMethod("genSeq", signature(obj = "bbcdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rbbcdSeq", 
                M = t(sapply(1:r, function(x) {
                  bbcdRand(N = N(obj), a = a(obj))
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
setMethod("genSeq", signature(obj = "bbcdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rbbcdSeq", 
                M = t(bbcdRand(N = N(obj), a = a(obj))), 
                a = a(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "bbcdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rbbcdSeq", 
                M = t(bbcdRand(N = N(obj), a = a(obj))), 
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
          signature(obj = "bbcdPar"),
          function(obj) {
            paste("bbcd(", obj@a, ")", sep = "")
          }
)

