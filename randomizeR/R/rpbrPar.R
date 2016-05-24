#' @include randPar.R
#' @include pbrPar.R
#' @include util.R
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
validaterpbrPar <- function(object) {
  errors <- character()
  rb <- randBlocks(object)
  ratio <- object@ratio
  filledBlock <- object@filledBlock
  
  if(!all(rb > 0)) {
    msg <- paste("At least one of the block lengths has value smaller or equal to zero. 
                 Should be greater than zero.")
    errors <- c(errors, msg)
  }
  
  if(!all(sapply(rb, function(x) x == round(x)))) {
    msg <- paste("At least one of the block lengths is not an integer. Should be though.")
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
# Class definition for pbrPar
# --------------------------------------------

# Randomization paramters class
setClass("rpbrPar",
         slots = c(rb = "numeric", filledBlock = "logical"),
         contains = "randPar",
         validity = validaterpbrPar)

#
# --------------------------------------------
# Constructor function for pbrPar
# --------------------------------------------

#' Representing Randomized Permuted Block Randomization
#' 
#' Represents the randomization procedure Randomized Permuted Block Randomization.
#'
#' @details
#' Fix the possible random block lengths \code{rb}, the number of treatment groups \code{K},
#' the sample size \code{N} and the vector of the \code{ratio}. Afterwards, one block length is
#' randomly selected of the random block lengths. The patients are assigned
#' according to the ratio to the corresponding treatment groups. This procedure is repeated
#' until \code{N} patients are assigned. Within each block all possbible
#' randomization sequences are equiprobable.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#' 
#' @return 
#' \code{S4} object of the class \code{rpbrPar}.
#'
#' @export
#'
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
rpbrPar <- function(N, rb, K = 2, ratio = rep(1, K), groups = LETTERS[1:K],
                    filledBlock = FALSE) {
  new("rpbrPar", rb = rb, N = N, K = K, ratio = ratio, groups = groups,
      filledBlock = filledBlock)
}


# --------------------------------------------
# Methods for pbrPar
# --------------------------------------------


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rpbrPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            bc <- genBlockConst(N(obj), randBlocks(obj), obj@filledBlock)
            new("rRpbrSeq", 
                M = t(blockRand(bc = bc, K = K(obj),
                ratio = ratio(obj))[1:N(obj)]), 
                filledBlock = obj@filledBlock,  
		            rb = randBlocks(obj),
                bc = list(bc),
                N = N(obj),
                K = K(obj),
		            ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rpbrPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            bc <- lapply(1:r, function(x) genBlockConst(N(obj),
                              randBlocks(obj), obj@filledBlock))
            new("rRpbrSeq", 
                M = t(sapply(bc, function(x) (blockRand(bc = x, K = K(obj),
                  ratio = ratio(obj)))[1:N(obj)])), 
                filledBlock = obj@filledBlock,  
		            rb = randBlocks(obj),
                bc = bc,				
                N = N(obj),
                K = K(obj),
		            ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rpbrPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            bc <- genBlockConst(N(obj), randBlocks(obj), obj@filledBlock)
            new("rRpbrSeq", 
                M = t(blockRand(bc = bc, K = K(obj), ratio = ratio(obj))[1:N(obj)]),
                filledBlock = obj@filledBlock,  
		            rb = randBlocks(obj),
		            bc = list(bc),
                N = N(obj),
                K = K(obj),
		            ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "rpbrPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            bc <- lapply(1:r, function(x) genBlockConst(N(obj),
                         randBlocks(obj), obj@filledBlock))
            new("rRpbrSeq", 
                M = t(sapply(bc, function(x) (blockRand(bc = x, K = K(obj),
                ratio = ratio(obj)))[1:N(obj)])), 
                filledBlock = obj@filledBlock,  
                bc = bc,
		            rb = randBlocks(obj),				
                N = N(obj),
                K = K(obj),
		            ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)


#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "rpbrPar"),
          function(obj) {
              rb <- capture.output(cat(obj@rb, sep = ","))
              if (obj@filledBlock) {
                paste(c("RPBRFB(", rb, ")"), sep = "", collapse = "")
              } else {
                paste(c("RPBR(", rb, ")"), sep = "", collapse = "")
              }  
          }
)
