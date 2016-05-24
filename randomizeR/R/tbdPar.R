#' @include randPar.R
#' @include util.R
NULL

###############################################
# --------------------------------------------#
# Class tbdPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the tbdPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatetbdPar <- function(object) {
  errors <- character()
  bc <- object@bc
  ratio <- object@ratio
  K <- object@K
  
  if(!all(bc %% sum(ratio) == 0)) {
    msg <- paste("One of the block length is not a multiple of sum(ratio) = "
                 , sum(ratio), ".", sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for tbdPar
# --------------------------------------------

# Randomization paramters generic 
setClass("tbdPar",
         slots = c(bc = "numeric"),
         contains = "randPar",
         validity = validatetbdPar)


# --------------------------------------------
# Constructor function for tbdPar
# --------------------------------------------

#' Representing Truncated Binomal Design
#' 
#' Represents the Truncated Binomial Design.
#'
#' @details
#' A fair toin is tossed until half of the patients have been assigned to one of
#' the treatment arms. Afterwards, the randomization list is filled with the
#' other treatment.
#'
#' @family randomization procedures
#' 
#' @inheritParams overview
#' 
#' @return
#' S4 object of the class \code{tbdPar}.
#' 
#' @export
#' 
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
tbdPar <- function(bc = N, groups = LETTERS[1:2]) {
  new("tbdPar", N = sum(bc), bc = bc, K = 2, ratio = c(1, 1), groups = groups)
}

# --------------------------------------------
# Sampling algorithm for TBD
# --------------------------------------------

# Truncated Binomial Design 
#
# This procedure generalises the Truncated Binomial Design.
#
# @inheritParams overview
# 
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
# is allocated to treatment A (resp. B).
# 
# @references 
# W. F. Rosenberger and J. M. Lachin: Randomization in Clinical Trials. Wiley
# (2002) 
tbdRand <- function(N, bc = N, K = 2, ratio = rep(1, K)) {
  stopifnot(all(is.numeric(bc)), all(bc %% 2 == 0))
  if(K > 2) stop("TBD: K>2 not available yet.")
  if(!all(ratio == 1)) stop("TBD: other ratios than 1:1 are not supported yet.") 
  ## generating the randomization sequence
  res <- sapply(bc, function(x) {
    R <- rbinom(x, 1, 0.5)
    lR <- lastRandom(R)
    if (R[lR] == 0) {
      R[(lR + 1):x] <- 1
    } else {
      R[(lR + 1):x] <- 0
    }
    
    R
  })
  as.vector(unlist(res))[1:N] 
}



# --------------------------------------------
# Methods for tbdPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "tbdPar"),
          function(obj) {
            if(obj@K != 2 || !identical(obj@ratio, c(1,1))) {
              stop("Only possible for K equals 2 and ratio corresponds to c(1,1).")
            }  
            allSeqs <- compltSet(obj)
            blockEnds <- cumsum(blocks(obj))
            bal <- apply(allSeqs,1, function(x, blockEnds) {
              all(cumsum(2*x - 1)[blockEnds] == 0)
            }, blockEnds = blockEnds)
            new("tbdSeq",
                M = allSeqs[bal, ],
                bc = blocks(obj),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "tbdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
	    set.seed(seed)
            new("rTbdSeq", 
                M = t(sapply(1:r,function(x) {
                  tbdRand(N(obj), blocks(obj), K(obj), ratio(obj))
                  })), 
                N = N(obj), 
                bc = obj@bc,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "tbdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rTbdSeq", 
                M = t(tbdRand(N(obj), blocks(obj), K(obj), ratio(obj))),
                N = N(obj),  
                bc = obj@bc,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "tbdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rTbdSeq", 
                M = t(sapply(1:r,function(x) {
                  tbdRand(N(obj), blocks(obj), K(obj), ratio(obj))
                  })), 
                N = N(obj), 
                bc = obj@bc,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "tbdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rTbdSeq", 
                M = t(tbdRand(N(obj), blocks(obj), K(obj), ratio(obj))),
                N = N(obj), 
                bc = obj@bc,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "tbdPar"),
          function(obj) {
            if (obj@N == obj@bc[1]) {
              "TBD"
            } else if (sum(!duplicated(obj@bc)) == 1) {
              paste("TBD(", obj@bc[1], ")", sep = "")
            } else {
              bc <- capture.output(cat(obj@bc, sep = ","))
              paste(c("TBD(", bc, ")"), sep = "", collapse = "")
            }  
          }
)
