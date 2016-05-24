#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class hadaPar                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for hadaPar
# --------------------------------------------

# Randomization paramters generic
setClass("hadaPar", contains = "randPar")

# --------------------------------------------
# Constructor function for hadaPar
# --------------------------------------------

#' Representing Hadamard Randomization
#' 
#' Represents the randomization procedure Hadamard Randomization.
#'
#' Hadamard randomization has been proposed by R.A. Bailey. The key idea is to
#' use the columns of a special Hadamard Matrix as a randomization scheme. The
#' implemented algorithm uses the Hadamard Matrix with \code{N=12} columns 
#' proposed in the paper, see references. 
#'
#' @note
#' \code{getProb} and \code{getAllSeq} are currently only supported for \code{hadaPar}
#' with total sample size \code{N=12}.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#' 
#' @return
#' \code{S4} object of the class \code{hadaPar}.
#'
#' @export
#'
#' @references 
#' R.A. Bailey and P.R. Nelson (2003) Hadamard Randomization: A valid restriction
#' of random permuted blocks. \emph{Biometrical Journal}, \strong{45}, 554-60.
hadaPar <- function(N, groups = LETTERS[1:2]) {
  new("hadaPar", N = N, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Methods for crSeq
# --------------------------------------------

# Hadamard Randomization
# 
# Computes a Hadamard Randomization sequence for a clinical trial with
# several blocks.
#
# @inheritParams overview
#
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
# is allocated to treatment A (resp. B).
# 
hadaRand <- function(bc) {
  stopifnot(is.numeric(bc))
  
  k <- c(1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,1,1,1,0,1,0,1,0,0,0,1,1,1,0,1,0,1,0,
         0,0,1,1,1,0,1,0,0,0,0,0,1,1,1,0,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,1,1,1,
         0,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,
         0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,0,0,1,0,0,0,1,1,1,0)
  H <- matrix(k, nrow = 11, ncol = 12)
  # nr of blocks needed
  nob <- ceiling(sum(bc)/12) 
  R <- numeric(0)
  for (i in 1:nob) {
    p <- H[sample(1:11, 1), ]
    # randomly select if 0 -> A (if 0) or 0 -> B (if 1)
    zeroToA <- sample(c(0, 1), 1) 
    p <- (zeroToA + p) %% 2
    R <- c(R, p)
  }
  R[1:sum(bc)]
}


# --------------------------------------------
# Methods for hadaPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "hadaPar"),
          function(obj) {
            if(obj@N > 12) stop("Only possible up to N equals 12.")
            k <- c(1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,1,1,1,0,1,0,1,0,0,0,1,1,1,
                   0,1,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,1,1,1,0,1,0,0,1,0,0,1,1,1,
                   0,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,1,1,1,0,1,0,0,1,0,0,0,1,1,
                   0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,
                   1,1,0,0,1,0,0,0,1,1,1,0)
            new("hadaSeq",
                M = matrix(k, nrow = 11, ncol = 12)[ ,1:N(obj)],
                N = N(obj),
                K = K(obj),
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "hadaPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rHadaSeq", 
                M = t(sapply(1:r, function(x) hadaRand(N(obj)))), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "hadaPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            if(K(obj) > 2) stop("HAD: K>2 not available.")
            new("rHadaSeq",
                M = t(hadaRand(N(obj))),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "hadaPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rHadaSeq", 
                M = t(sapply(1:r, function(x) hadaRand(N(obj)))), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "hadaPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            if(K(obj) > 2) stop("HAD: K>2 not available.")
            new("rHadaSeq",
                M = t(hadaRand(N(obj))),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "hadaPar"),
          function(obj) {
            "HADA"
          }
)