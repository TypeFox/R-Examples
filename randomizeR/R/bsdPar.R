#' @include randPar.R
#' @include ebcPar.R
NULL

###############################################
# --------------------------------------------#
# Class bsdPar                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the bsdPar class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validatebsdPar <- function(object) {
  errors <- character()
  mti <- object@mti
  lengthMTI <- length(mti)

  if(lengthMTI != 1) {
    msg <- paste("mti has length ", lengthMTI, ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if (round(mti[1]) != mti) {
    msg <- paste("First element of mti is ", mti, ". Should be an integer.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(mti[1] < 0){
    msg <- "mti must be a positive integer"
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for bsdPar
# --------------------------------------------

# Randomization paramters generic
setClass("bsdPar",
         slots = c(mti = "numeric"),
         contains = "randPar",
         validity = validatebsdPar)

# --------------------------------------------
# Constructor function for bsdPar
# --------------------------------------------

#' Representing Big Stick Design
#' 
#' Represents the randomization procedure Big Stick Design.
#'
#' @details
#' Tossing a fair coin as long as the difference in group sizes doesn`t
#' exceed the \code{mti}. If the \code{mti} is reached a deterministic
#' allocation is done, so that the difference in group sizes is reduced.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @references
#' J. F. Soares and C. F. Jeff Wu (1983) Some Restricted Randomization Rules in
#' Sequential Designs. \emph{Comm. in Stat.}, \strong{12}, 2017-34. 
#' 
#' @return 
#' \code{S4} object of the class \code{bsdPar}.
#'
#' @export
bsdPar <- function(N, mti, groups = LETTERS[1:2]) {
  new("bsdPar", N = N, mti = mti, K = 2, ratio = c(1, 1), groups = groups)
}


# --------------------------------------------
# Sampling algorithm for EBC
# --------------------------------------------

#' Sampling algorithm for BSD 
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
#' J. F. Soares and C. F. Jeff Wu (1983) Some Restricted Randomization Rules in
#' Sequential Designs. Comm. in Stat., 12, 2017-34. 
bsdRand <- function(N, mti, K = 2) {
  stopifnot(round(N) == N, round(mti) == mti)
  if (!(K == 2)) stop("EBC, BSD: K>2 not available yet.")
  R <- numeric(N); reps <- 0; sumR <- 0; imb <- 0
  while(reps < N) {
    # case analysis
    if (abs(imb) < mti) {
      R[reps + 1] <- rbinom(1, 1, 0.5)  
    } else if (imb == mti) {
      R[reps + 1] <- 0
    }  else {
      R[reps+1] <- 1  
    }
    reps <- reps + 1
    sumR <- sumR + R[reps]
    imb <- 2*sumR - reps
  }
  R
}



# --------------------------------------------
# Methods for bsdPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "bsdPar"),
          function(obj) {
            allSeqs <- compltSet(obj)
            inside <- apply(allSeqs,1, function(x, mti) {
                            all(abs(cumsum(2*x - 1)) <= mti)
                            }, mti = mti(obj))
            new("bsdSeq", 
		          M = allSeqs[inside, ], 
		          mti = mti(obj), 
		          N = N(obj),
             		  K = K(obj),
		          ratio = obj@ratio,
              		  groups = obj@groups
		          )
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "bsdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rBsdSeq", 
                M = t(sapply(1:r, function(x) {
                  bsdRand(N = N(obj), mti = mti(obj))
                  })),
                N = N(obj),
                mti = mti(obj),
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)
#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "bsdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rBsdSeq", 
              M = t(sapply(1:r, function(x) {
                  bsdRand(N = N(obj), mti = mti(obj))
                  })),
              N = N(obj),
              mti = mti(obj),
              K = K(obj),
              groups = obj@groups,
              ratio = obj@ratio,
	            seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "bsdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rBsdSeq", 
                M = t(bsdRand(N = N(obj), mti = mti(obj))), 
                mti = mti(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
		            seed = seed)
          }
)


#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "bsdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rBsdSeq", 
                M = t(bsdRand(N = N(obj), mti = mti(obj))), 
                mti = mti(obj), 
                N = N(obj),
                K = K(obj),
                groups = obj@groups,
                ratio = obj@ratio,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "bsdPar"),
          function(obj) {
            paste("BSD(", obj@mti, ")", sep = "")
          }
)

