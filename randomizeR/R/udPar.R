#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class udPar                                 #
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
validateudpar <- function(object) {
  errors <- character()
  ini <- object@ini
  add <- object@add

  if(length(ini) != 1) {
    msg <- paste("ini has length  ", length(ini), ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if (round(ini[1]) != ini) {
    msg <- paste("First element of ini is  ", ini, ". Should be an integer.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if (ini[1] < 0) {
    msg <- paste("First element of ini is negative should be positive.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(length(add) != 1) {
    msg <- paste("add has length  ", length(add), ". Should be length one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if (round(add[1]) != add) {
    msg <- paste("First element of add is  ", add, ". Should be an integer.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if (add[1] <= 0) {
    msg <- paste("First element of add is negative or zero should be positive.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for udPar
# --------------------------------------------

# Randomization paramters generic 
setClass("udPar",
         slots = c(ini = "numeric", add = "numeric"),
         contains = "randPar",
         validity = validateudpar)


# --------------------------------------------
# Constructor function for udPar
# --------------------------------------------

#' Representing Wei's Urn Design
#' 
#' Represents Wei's Urn Design.
#'
#' @details
#' An urn is filled with a number of \code{ini} balls of both of the treatments.
#' Afterwards, a ball is drawn randomly from the urn. Finally, \code{add} balls
#' are added to the urn from the opposite treatment. This procedure is repeated until
#' \code{N} patients are assigend.
#' 
#' @family randomization procedures
#' 
#' @inheritParams overview
#'
#' @return 
#' \code{S4} object of the class \code{udPar}.
#'
#' @export
#'
#' @references
#' L.J. Wei (1977) A Class of Designs for Sequential Clinical Trials.
#' \emph{Journal of the American Statistical Association}, \strong{72}, 382-6.
udPar <- function(N, ini, add, groups = LETTERS[1:2]) {
  new("udPar", N = N, ini = ini, add = add, K = 2, ratio = c(1, 1),
      groups = groups)
}


# --------------------------------------------
# Sampling algorithm for UD
# --------------------------------------------

# Wei's Urn design
# 
# Computes a randomisation sequence based on Weis Urn Design
#
# @inheritParams overview
# @param ini integer representing the initial urn composition.
# @param add integer representing the number of balls that are added to the
# urn in each step.
#
# @return A vector with the allocation sequence for a clinical trial. 
# It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
# is allocated to treatment A (resp. B).
# 
# @export
urnRand <- function(N, ini, add) {
  R <- numeric(N)
  R[1] <- sample(c(0, 1), size = 1)
  sumR <- R[1]
  for(j in 1:(N-1)) {
    R[j + 1] <- rbinom(1, 1, prob = (ini + add*(j-sumR))/(2*ini + add*j))
    sumR <- sumR + R[j+1]
    }
  R
}


# --------------------------------------------
# Methods for udPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "udPar"),
          function(obj) {
                if(obj@ini == 0) {
                  allSeqs <- compltSet(obj)
                  inside <- apply(allSeqs, 1, function(x) {
                    if(sum(x[1:2]) %in% c(0, 2)) {
                       FALSE
                    } else {
                      TRUE
                    }
                  })
                  M = allSeqs[inside, ]
                } else {
                  M = compltSet(obj)
                }
                new("udSeq",
                    M = M, 
                    N = N(obj),
                    ini = obj@ini, 
                    add = obj@add,
                    ratio = obj@ratio,
                    groups = obj@groups, 
                    K = K(obj))
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "udPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
	    set.seed(seed)
            new("rUdSeq", 
                M = t(sapply(1:r, function(x) {
                  urnRand(N(obj), obj@ini, obj@add)
                  })), 
                N = N(obj),
                ini = obj@ini, 
                add = obj@add,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "udPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            if(K(obj) > 2) stop("UD: K>2 not available.")
            new("rUdSeq", 
                M = t(urnRand(N(obj), obj@ini, obj@add)), 
                N = N(obj), 
                ini = obj@ini, 
                add = obj@add,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)


#' @rdname generateRandomSequences 
setMethod("genSeq", signature(obj = "udPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
	    seed <- sample(.Machine$integer.max, 1)
	    set.seed(seed)
            new("rUdSeq", 
                M = t(sapply(1:r,function(x) {
                  urnRand(N(obj), obj@ini, obj@add)
                  })), 
                N = N(obj),
                ini = obj@ini, 
                add = obj@add,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "udPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            if(K(obj) > 2) stop("UD: K>2 not available.")
            new("rUdSeq", 
                M = t(urnRand(N(obj), obj@ini, obj@add)), 
                N = N(obj), 
                ini = obj@ini, 
                add = obj@add,
                K = K(obj),
                ratio = obj@ratio,
                groups = obj@groups,
		            seed = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "udPar"),
          function(obj) {
            paste("UD(", obj@ini, ",", obj@add, ")", sep = "")
          }
)

