#' @include randPar.R
#' @include tbdPar.R
#' @include util.R
NULL

###############################################
# --------------------------------------------#
# Class rTbdPar                                #
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
validatertbdPar <- function(object) {
  errors <- character()
  rb <- randBlocks(object)
  ratio <- object@ratio
  K <- object@K
  filledBlock <- object@filledBlock
  
  if(!all(rb > 0)) {
    msg <- paste("At least one of the block lengths has value smaller or equal to zero.",  
                 "Should be greater than zero.")
    errors <- c(errors, msg)
  }
  
  if(!all(rb %% K == 0)) {
    msg <- paste("All block lengths should be even integers.",  
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(!all(rb %% sum(ratio) == 0)) {
    msg <- paste("One of the block length is not a multiple of sum(ratio) = "
                 , sum(ratio), ".", sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(filledBlock) > 1) {
    msg <- paste("filledBlock has length  ", length(filledBlock), ". Should be one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for tbdPar
# --------------------------------------------

# Randomization paramters generic 
setClass("rtbdPar",
         slots = c(rb = "numeric", filledBlock = "logical"),
         contains = "randPar",
         validity = validatertbdPar)


# --------------------------------------------
# Constructor function for tbdPar
# --------------------------------------------

#' Representing Randomized Truncated Binomial Design
#' 
#' Represents the randomization procedure Randomized Truncated Binomial Design.
#'
#' @details
#' Fix the possible random block lengths \code{rb} and the sample size of the trial \code{N}.
#'  Afterwards, one block length is
#' randomly selected of the random block lengths. In this block a fair coin is tossed
#' for the patient assignments until half of the patients have been assigend to one of
#' the treatment arms. Afterwards, the block is filled with the
#' other treatment. This procedure is repeated
#' until \code{N} patients are assigned.
#'
#' @family randomization procedures
#' 
#' @inheritParams overview
#' 
#' @return
#' \code{S4} object of the class \code{rtbdPar}.
#' 
#' @export
#' 
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
rtbdPar <- function(N, rb = N, groups = LETTERS[1:2], filledBlock = FALSE){
  new("rtbdPar", N = N, rb = rb, K = 2, ratio = c(1, 1), groups = groups,
      filledBlock = filledBlock)
}


# --------------------------------------------
# Methods for rtbdPar
# --------------------------------------------

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rtbdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
	    set.seed(seed)
	    bc <- lapply(1:r, function(x) genBlockConst(N(obj),
              randBlocks(obj), obj@filledBlock))
            new("rRtbdSeq", 
               M = t(sapply(bc, function(x) tbdRand(N(obj), x, K(obj),
                 ratio(obj)))), 
               filledBlock = obj@filledBlock, 
               N = N(obj), 
               rb = randBlocks(obj),
      	       bc = bc,
               K = K(obj),
               ratio = obj@ratio,
               groups = obj@groups,
	             seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rtbdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
	    set.seed(seed)
	    bc <- genBlockConst(N(obj), randBlocks(obj), obj@filledBlock)
      new("rRtbdSeq", 
         M = t(tbdRand(N(obj), bc, K(obj), ratio(obj))),
         filledBlock = obj@filledBlock, 
         N = N(obj), 
         rb = randBlocks(obj),
         bc = list(bc),
         K = K(obj),
         ratio = obj@ratio,
         groups = obj@groups,
         seed = seed)
      }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rtbdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
		seed <- sample(.Machine$integer.max, 1)
		set.seed(seed)
		bc <- lapply(1:r, function(x) genBlockConst(N(obj),
                  randBlocks(obj), obj@filledBlock))
                new("rRtbdSeq", 
                  M = t(sapply(bc, function(x) tbdRand(N(obj), x, K(obj),
                    ratio(obj)))), 
                  filledBlock = obj@filledBlock, 
                  N = N(obj), 
                  rb = randBlocks(obj),
		              bc = bc,
                  K = K(obj),
		              ratio = obj@ratio,
                  groups = obj@groups,
		              seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rtbdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
	    bc <- genBlockConst(N(obj), randBlocks(obj), obj@filledBlock)
      new("rRtbdSeq", 
          M = t(tbdRand(N(obj), bc, K(obj), ratio(obj))),
          filledBlock = obj@filledBlock, 
          N = N(obj),  
          rb = randBlocks(obj),
          bc = list(bc),
          K = K(obj),
          ratio = obj@ratio,
          groups = obj@groups,
          seed = seed)
     }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "rtbdPar"),
          function(obj) {
            rb <- capture.output(cat(obj@rb, sep = ","))
            if (obj@filledBlock) {
              paste(c("RTBDFB(", rb, ")"), sep = "", collapse = "")
            } else {
              paste(c("RTBD(", rb, ")"), sep = "", collapse = "")
            }
          }
)
