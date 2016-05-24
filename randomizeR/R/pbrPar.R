#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class pbrPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the pbrPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatepbrPar <- function(object) {
  errors <- character()
  bc <- object@bc
  ratio <- object@ratio
  
#   if(!all(bc > 0)) {
#     msg <- paste("At least one of the block lengths has value smaller or equal to zero. 
#                  Should be greater than zero.")
#     errors <- c(errors, msg)
#   }
#   
#   if(!all(sapply(bc, function(x) x == round(x)))) {
#     msg <- paste("At least one of the block lengths is not an integer. Should be though.")
#     errors <- c(errors, msg)
#   }
  
  if(!all(bc %% sum(ratio) == 0)) {
    msg <- paste("One of the block lengths is not a multiple of sum(ratio) = "
                  , sum(ratio), ".", sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for pbrPar
# --------------------------------------------

# Randomization paramters generic
setClass("pbrPar",
         slots = c(bc = "numeric"),
         contains = "randPar",
         validity = validatepbrPar)


# --------------------------------------------
# Constructor function for pbrPar
# --------------------------------------------

#' Representing Permuted Block Randomization
#' 
#' Represents the randomization procedure Permuted Block Randomization.
#'
#' @details
#' Fix the block constellation \code{bc}, the number of treatment groups \code{K},
#' and the vector of the \code{ratio}. Afterwards,
#' in each block the patients are assigned according to the ratio to the 
#' corresponding treatment groups. All generated randomization sequences
#' are equiprobable.
#' 
#' @family randomization procedures 
#' 
#' @inheritParams overview
#' 
#' @return 
#' \code{S4} object of the class \code{pbrPar}.
#'
#' @export
#'
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
pbrPar <- function(bc, K = 2, ratio = rep(1, K), groups = LETTERS[1:K]) {
  new("pbrPar", bc = bc, N = sum(bc), K = K, ratio = ratio, groups = groups)
}


# --------------------------------------------
# Sampling alogorithm for PBR
# --------------------------------------------

# Permuted block randomization
#
# Compute a permuted block randomization sequence for a clinical trial with
# several blocks.
#
# @inheritParams overview
# 
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
# is allocated to treatment A (resp. B).
# 
blockRand <- function(bc, K = 2, ratio = rep(1, K)) {
  as.vector(unlist(sapply(bc, function(k) blockSeq(k, K, ratio))))
}

# Permuted block randomization
# 
# Compute a permuted block randomization sequence for a clinical trial for
# one block.
#
# @inheritParams overview
# 
# @return 
# A vector with the allocation sequence for a clinical trial for one block.
# 
blockSeq <- function(k, K = 2, ratio = rep(1, K)) {
  ratio <- ratio/sum(ratio)
  sample(rep(0:(K-1), times = ratio*k))
}


# --------------------------------------------
# Methods for pbrPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", 
          signature(obj = "pbrPar"),
          function(obj) {
            if(obj@K != 2 || !identical(obj@ratio, c(1,1))) {
              stop("Only possible for K equals 2 and ratio corresponds to c(1,1).")
            }  
            allSeqs <- compltSet(obj)
            blockEnds <- cumsum(blocks(obj))
            bal <- apply(allSeqs,1, function(x, blockEnds) {
              all(cumsum(2*x-1)[blockEnds] == 0)
              }, blockEnds = blockEnds)
            new("pbrSeq",
                M = allSeqs[bal, ],
                bc = blocks(obj),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "pbrPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rPbrSeq", 
                M = t(blockRand(bc = blocks(obj), K = K(obj),
                  ratio = ratio(obj))), 
                bc = blocks(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
                seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "pbrPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rPbrSeq", 
                M = t(sapply(1:r, function(x) {
                  blockRand(bc = blocks(obj), K = K(obj), ratio = ratio(obj))
                  })), 
                bc = blocks(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "pbrPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rPbrSeq", 
                M = t(blockRand(bc = blocks(obj), K = K(obj),
                ratio = ratio(obj))), 
                bc = blocks(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", 
          signature(obj = "pbrPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rPbrSeq", 
                M = t(sapply(1:r, function(x) {
                  blockRand(bc = blocks(obj), K = K(obj), ratio = ratio(obj))
                  })), 
                bc = blocks(obj), 
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "pbrPar"),
          function(obj) {
            if (sum(!duplicated(obj@bc)) == 1) {
              paste("PBR(", obj@bc[1], ")", sep = "")
            } else {
              bc <- capture.output(cat(obj@bc, sep = ","))
              paste(c("PBR(", bc, ")"), sep = "", collapse = "")
            }  
          }
)


