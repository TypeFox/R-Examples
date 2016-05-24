#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class gbcdPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the gbcdPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validategbcdPar <- function(object) {
  errors <- character()
  rho <- object@rho
  lengthRHO <- length(rho)
  
  if(lengthRHO != 1) {
    msg <- paste("rho has length ", lengthRHO, ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if (rho[1] < 0) {
    msg <- paste("First element of rho is ", rho, ". Should be positive.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for gbcdPar
# --------------------------------------------

# Randomization paramters generic
setClass("gbcdPar",
         slots = c(rho = "numeric"),
         contains = "randPar",
         validity = validategbcdPar)

# --------------------------------------------
# Constructor function for gbcdPar
# --------------------------------------------

#' Representing Generalized Biased Coin Design
#' 
#' Represents the randomization procedure Generalized Biased Coin Design.
#'
#' @details
#' Generalization of Wei's urn and Efron's biased coin design.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @references
#' R. L. Smith (1984) Sequential treatment allocation using biased coin designs. 
#' \emph{Journal of the Royal Statistical Society B},
#' \strong{46}, 519-543. \cr
#' W. F. Rosenberger and J. M. Lachin (2002) Randomization in Clinical Trials. \emph{Wiley},
#' 64-65
#' 
#' @return 
#' \code{S4} object of the class \code{gbcdPar}.
#'
#' @export
gbcdPar <- function(N, rho, groups = LETTERS[1:2]) {
  new("gbcdPar", N = N, rho = rho, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Sampling algorithm for GBCD
# --------------------------------------------

#' Sampling algorithm for gbcd 
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
#' R. L. Smith (1984) Sequential treatment allocation using biased coin designs. 
#' \emph{Journal of the Royal Statistical Society B},
#' \strong{46}, 519-543. \cr
#' W. F. Rosenberger and J. M. Lachin (2002) Randomization in Clinical Trials. \emph{Wiley},
#' 64-65
gbcdRand <- function(N, rho, K = 2) {
  stopifnot(round(N) == N)
  if (!(K == 2)) stop("GBCD: K > 2 not available yet.")
  R <- numeric(N); reps <- 0; nA <- 0; nB <- 0; d <- 0
  while(reps < N) {
    # case analysis
    if (abs(d) == 0) {
      R[reps + 1] <- rbinom(1, 1, 0.5)  
    } else {
      R[reps + 1] <- rbinom(1, 1, 1 - phi)  
    }
    reps <- reps + 1
    nB <- sum(R[1:reps])
    nA <- reps - nB
    d <- nA-nB
    phi <- nB^rho/(nA^rho + nB^rho)
  }
  R
}



# --------------------------------------------
# Methods for gbcdPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "gbcdPar"),
          function(obj) {
            allSeqs <- compltSet(obj) 
            # Two successive treatment assignments at the beginning of the trial 
            # cannot be the same 
            inside <- which(apply(allSeqs, 1, function(x) sum(x[1:2])==1)==TRUE)
            new("gbcdSeq", 
                M = allSeqs[inside,],
                rho = rho(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups
            )
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "gbcdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rGbcdSeq", 
                M = t(sapply(1:r, function(x) {
                  gbcdRand(N = N(obj), rho = rho(obj))
                })),
                N = N(obj),
                rho = rho(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
                seed = seed)
          }
)
#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "gbcdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rGbcdSeq", 
                M = t(sapply(1:r, function(x) {
                  gbcdRand(N = N(obj), rho = rho(obj))
                })),
                N = N(obj),
                rho = rho(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "gbcdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rGbcdSeq", 
                M = t(gbcdRand(N = N(obj), rho = rho(obj))), 
                rho = rho(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "gbcdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rGbcdSeq", 
                M = t(gbcdRand(N = N(obj), rho = rho(obj))), 
                rho = rho(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "gbcdPar"),
          function(obj) {
            paste("gbcd(", obj@rho, ")", sep = "")
          }
)

