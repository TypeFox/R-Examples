#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class chenPar                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the chenPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateChenPar <- function(object) {
  errors <- character()
  mti <- object@mti
  p <- object@p
  lengthMTI <- length(mti)
  
  if(lengthMTI != 1) {
    msg <- paste("mti has length ", lengthMTI, ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(round(mti[1]) != mti) {
    msg <- paste("First element of mti is ", mti, ". Should be an integer.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(mti[1] < 0){
    msg <- "mti must be a positive integer"
    errors <- c(errors, msg)
  }
  
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
# Class definition for chenPar
# --------------------------------------------

# Randomization paramters generic
setClass("chenPar",
         slots = c(mti="numeric", p = "numeric"),
         contains = "randPar",
         validity = validateChenPar)

# --------------------------------------------
# Constructor function for chenPar
# --------------------------------------------

#' Representing Chen's Design
#' 
#' Represents the randomization procedure Chen's Design.
#'
#' @details
#' Flip a biased coin with probability \code{p} in favour of the treatment
#' which is allocated less frequently as long as the difference in group sizes doesn`t
#' exceed the \code{mti}. If the \code{mti} is reached a deterministic
#' allocation is done, so that the difference in group sizes is reduced.
#' If both treatments have been assigned equally often a fair coin is tossed.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @references
#' Chen Yung-Pin (1999) Biased coin design with imbalance tolerance.
#' \emph{Comm. in Stat.}, \strong{15}, 953-975. 
#' 
#' @return 
#' \code{S4} object of the class \code{chenPar}.
#'
#' @export
chenPar <- function(N, mti = N, p = 0.5, groups = LETTERS[1:2]) {
  new("chenPar", N = N, mti = mti, p = p, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Sampling algorithm for Chen
# --------------------------------------------

#' Representing Chen's Design
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
#' Chen Yung-Pin (1999) Biased coin design with imbalance tolerance.
#' Comm. in Stat., 15, 953-975. 
chenRand <- function(N, mti, p, K = 2) {
  stopifnot(round(N) == N, round(mti) == mti, is.numeric(p), p >= 0.5 & p <= 1)
  if (!(K == 2)) stop("EBC, BSD, CHEN: K>2 not available yet.")
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
# Methods for chenPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "chenPar"),
          function(obj) {
            allSeqs <- compltSet(obj)
            inside <- apply(allSeqs,1, function(x, mti) {
              all(abs(cumsum(2*x - 1)) <= mti)
            }, mti = mti(obj))
            new("chenSeq", 
                M = allSeqs[inside, ], 
                mti = mti(obj),
                p = coin(obj),
                N = N(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups
            )
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "chenPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rChenSeq", 
                M = t(sapply(1:r, function(x) {
                  chenRand(N = N(obj), mti = mti(obj), p = coin(obj))
                })),
                N = N(obj),
                mti = mti(obj),
                p = coin(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
                seed = seed)
          }
)
#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "chenPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rChenSeq", 
                M = t(sapply(1:r, function(x) {
                  chenRand(N = N(obj), mti = mti(obj), p = coin(obj))
                })),
                N = N(obj),
                mti = mti(obj),
                p = coin(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "chenPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rChenSeq", 
                M = t(chenRand(N = N(obj), mti = mti(obj), p = coin(obj))), 
                mti = mti(obj), 
                p = coin(obj),
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "chenPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rChenSeq", 
                M = t(chenRand(N = N(obj), mti = mti(obj), p = coin(obj))), 
                mti = mti(obj), 
                p = coin(obj),
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
                seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "chenPar"),
          function(obj) {
            paste("CHEN(", obj@mti, ", " ,  round(obj@p, digits = 2), ")", sep = "")
          }
)

