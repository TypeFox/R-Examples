#' @include getDesign.R
NULL

###############################################
# --------------------------------------------#
# Class randPar                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validate randomization generators
#
# @param object object 
validateRandPar <- function(object) {
  errors <- character()
  N <- object@N
  K <- object@K
  groups <- object@groups
  ratio <- object@ratio

  if(N <= 1) {
    msg <- paste("N should be greater than one.",   
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(!(length(N) == 1)) {
    msg <- paste("N = ", N, " should have length 1. Has length ", length(N),
               ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(N == ceiling(N))) {
    msg <- paste("N should be an integer.", 
                  sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(length(K) == 1)) {
    msg <- paste("K =", K, " should have length 1. Has length ", length(K),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(K == ceiling(K))) {
    msg <- paste("K should be an integer.", 
                  sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(length(groups) == K) && K == ceiling(K)) {
    msg <- paste("Length of groups is ", length(groups), ". Should have length ", K,
                 "." , sep = "")
    errors <- c(errors, msg)
  }

  if(sum(duplicated(groups)) > 0) {
    msg <- paste("Duplicated group names selected, must be unique.",
                 sep = "")
    errors <- c(errors, msg)
  }

  if(!(length(ratio) == K) && K == ceiling(K)) {
    msg <- paste("Length of ratio is ", length(ratio), ". Should have length ", K,
                 "." , sep = "")
    errors <- c(errors, msg)
  }

  if(!(all(ratio >= 1))) {
    msg <- paste("All entries of ratio must be greater than 1.",
                 sep = "")
    errors <- c(errors, msg)
  }

  if(!(all(ratio == ceiling(ratio)))) {
    msg <- paste("All entries of ratio must be integers.",
                 sep = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for randPar
# --------------------------------------------

#' Randomization paramters generic
setClass("randPar", 
         slots = c(N = "numeric", K = "numeric" , ratio = "numeric",
           groups = "character"),
         validity = validateRandPar)

# --------------------------------------------
# Constructor function for randPar
# --------------------------------------------

#' Settings for randomization procedures
#' 
#' Randomization procedures in randomizeR are represented by objects that inherit 
#' from \code{randPar}. The representation can then be used in order to
#' generate randomization sequences. In order generate a representation of a 
#' randomization procedure, call \code{\link{createParam}} or one of the following
#' functions.
#' 
#' @section Supported randomization procedures: 
#' \itemize{
#'   \item Complete Randomization (\code{\link{crPar}})
#'   \item Efron's Biased Coin Design (\code{\link{ebcPar}})
#'   \item Generalized Biased Coin Design (\code{\link{gbcdPar}})
#'   \item Adjustable Biased Coin Design (\code{\link{abcdPar}})
#'   \item Bayesian Biased Coin Design (\code{\link{bbcdPar}})
#'   \item Hadamard Randomization (\code{\link{hadaPar}})
#'   \item Maximal Procedure (\code{\link{mpPar}})
#'   \item Permuted Block Randomization (\code{\link{pbrPar}})
#'   \item Random Allocation Rule (\code{\link{rarPar}})
#'   \item Permuted Block Randomization with random block length (\code{\link{rpbrPar}})
#'   \item Truncated Binomial Design with random block length (\code{\link{rtbdPar}})
#'   \item Truncated Binomial Design (\code{\link{tbdPar}})
#'   \item Wei's Urn Design (\code{\link{udPar}})
#'   \item Chen's Design (\code{\link{chenPar}})
#' }
#' 
#' @seealso Generate randomization sequences \code{\link{genSeq}}.
#' Calculate the the complete set of randomization sequences of a randomization 
#' procedure.
#' \code{\link{getAllSeq}}.
#' 
#' @inheritParams overview
#'
#' @name randPar
NULL


# --------------------------------------------
# Accesssor functions for randPar
# --------------------------------------------

#' Method defining the $ operator for the randPar class
#' 
#' @inheritParams overview
setMethod("$", "randPar",
          function(x, name) slot(x, name))

#' Function returning the sample size slot of an S4 object
#'
#' @param obj object inheriting from randPar 
#' 
#' @export
N <- function(obj) {
  if (.hasSlot(obj, "N")) {
    obj@N
  } else {
    stop("object has no slot named N.")
  }
}

#' Function returning the block slot of an S4 object
#'
#' @param obj object of class pbrPAr
#' 
#' @export
blocks <- function(obj) {
  if (.hasSlot(obj, "bc")) {
    return(obj@bc)
  } else {
    stop("object has no slot for blocks.")
  }
}

#' Function returning the block slot of an S4 object
#'
#' @param obj object of class pbrPAr
#' 
#' @export
randBlocks <- function(obj) {
  if (.hasSlot(obj, "rb")) {
    obj@rb
  } else {
    stop("object has no random blocks.")
  }
}
  
#' Function returning the MTI slot of an S4 object
#'
#' @param obj object of class bsdPar or mpPar
#' 
#' @export
mti <- function(obj) {
  if (.hasSlot(obj, "mti")) {
    obj@mti
  } else {
    stop("object has no slot named mti.")
  }
}

#' Function returning the coin slot of an S4 object
#'
#' @param obj object extending class randPar or randSeq
#' @export
coin <- function(obj) {
  if (.hasSlot(obj, "p")) {
    obj@p
  } else stop("object has no slot named p.") 
}

#' Function returning the total sample size slot of an S4 object
#'
#' @param obj object of class randPar
#' 
#' @export
K <- function(obj) {
  if (.hasSlot(obj, "K")) {
    obj@K
  } else {
    stop("object has no slot named K.")
  }
}

#' Function returning the allocation ratio slot of an S4 object
#'
#' @param obj object of class randPar 
#' 
#' @export
ratio <- function(obj) {
  if (.hasSlot(obj, "ratio")) {
    obj@ratio
  } else {
    stop("object has no slot named ratio.")
  }
}

#' Function returning the allocation ratio slot of an S4 object
#'
#' @param obj object of class randPar
#' 
#' @export
method <- function(obj) {
  toupper(sub("Par", "", class(obj)[1]))
}

#' Function returning the adjusting parameter rho slot of an S4 object
#'
#' @param obj object of class randPar 
#' 
#' @export
rho <- function(obj) {
  if (.hasSlot(obj, "rho")) {
    obj@rho
  } else {
    stop("object has no slot named rho.")
  }
}

#' Function returning the adjusting parameter a slot of an S4 object
#'
#' @param obj object of class randPar 
#' 
#' @export
a <- function(obj) {
  if (.hasSlot(obj, "a")) {
    obj@a
  } else {
    stop("object has no slot named a.")
  }
}

# --------------------------------------------
# Show function for randPar
# --------------------------------------------

setMethod("show", "randPar", function(object) {
  validObject(object)
  # headline
  cat("\nObject of class \"", class(object)[1], "\"\n\n", sep = "")
  # crop the method from the class name of the randPar object
  cat("design =", getDesign(object), "\n") 
  # iterate through all slots of the randPar object
  names <- slotNames(object) 
  if (K(object) == 2) names <- names[!(names %in% "K")]
  if (all(ratio(object) == rep(1, K(object)))) {
    names <- names[!(names %in% "ratio")]
  }
  
  for(name in names) {
    cat(name, "=", slot(object, name),"\n")
  }
  cat("\n") 
})


# --------------------------------------------
# Generic functions for randPar
# --------------------------------------------

#' Complete set of randomization sequences
#' 
#' Outputs all randomization sequences for the given randomization procedure 
#' along with the parameters belonging to the randomization procedure.
#' The output consists of the parameters used for the generation of the 
#' randomization sequences (see \code{\link{createParam}}) and the matrix \code{M}
#' that stores the randomization sequences in its rows.
#' 
#' @details \code{getAllSeq} is a generic function which dispatches different 
#' methods depending on the type of input. 
#' 
#' @inheritParams overview
#' 
#' @return An object inheriting from \linkS4class{randSeq}, representing the set 
#' of randomization sequences for the given parameters.
#' The output consists of the parameters used for the generation of the 
#' randomization sequences (see \code{\link{createParam}}) and the matrix \code{M}
#' that stores the randomization sequences in its rows.
#' 
#' @seealso \code{\link{createParam}}
#' 
#' @examples
#' # CR
#' myPar <- crPar(6)
#' getAllSeq(myPar)
#' 
#' # EBC
#' myPar <- ebcPar(6, 0.667)
#' getAllSeq(myPar)
#' 
#' # BSD
#' myPar <- bsdPar(6, 2)
#' getAllSeq(myPar)
#' 
#' # PBR
#' myPar <- pbrPar(c(4, 2))
#' getAllSeq(myPar)
#' 
#' # RAR
#' myPar <- rarPar(8)
#' getAllSeq(myPar)
#' 
#' # MP 
#' myPar <- mpPar(8, 2)
#' getAllSeq(myPar)
#' 
#' # HAD
#' myPar <- hadaPar(8)
#' getAllSeq(myPar)
#' 
#' # TBD
#' myPar <- tbdPar(8)
#' getAllSeq(myPar)
#' 
#' # GBCD
#' myPar <- gbcdPar(8, 2)
#' getAllSeq(myPar)
#' 
#' # ABCD
#' myPar <- abcdPar(8, 3)
#' getAllSeq(myPar)
#'
#' # BBCD
#' myPar <- bbcdPar(8, 5)
#' getAllSeq(myPar)
#' 
#' # CHEN
#' myPar <- chenPar(8, 2, 0.667)
#' getAllSeq(myPar)
#' 
#' @name generateAllSequences
NULL

#' @rdname generateAllSequences
#'
#' @export
setGeneric("getAllSeq", function(obj) standardGeneric("getAllSeq"))

#' Generate random sequences
#' 
#' Generates a randomization sequences for a given randomization procedure.
#' 
#' @details
#' \code{genSeq} generates randomization sequences for a randomization 
#' procedure as defined by the input paramters. \code{genSeq} has two modes, 
#' according to the input.
#' \enumerate{
#'   \item \code{genSeq(obj,r)}: gives \code{r} random sequences from the 
#'   design specified by \code{obj}, along with the parameters stored in \code{obj}.
#'   \item \code{genSeq(obj)}: gives one random sequences from the 
#'   design specified by \code{obj}, along with the parameters stored in \code{obj}.
#' }
#' 
#' @inheritParams overview
#' 
#' @return An object inheriting from \linkS4class{randSeq}, representing the \code{r}
#' randomisation sequences generated at random for the specified randomization procedure.
#' The output consists of the parameters used for the generation of the 
#' randomization sequences (see \code{\link{createParam}}) and the matrix \code{M}
#' that stores the randomization sequences in its \code{r} rows.
#' If \code{r} is missing, one sequence is generated by default.
#' 
#' @examples
#' # CR
#' myPar <- crPar(10)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # EBC
#' myPar <- ebcPar(10, 0.667)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # BSD
#' myPar <- bsdPar(10, 2)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # PBR
#' myPar <- pbrPar(c(4, 4))
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # RAR
#' myPar <- rarPar(10)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # MP 
#' myPar <- mpPar(10, 2)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # HAD
#' myPar <- hadaPar(10)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # UD
#' myPar <- udPar(8, 0, 1)
#' genSeq(myPar,4)
#' genSeq(myPar)
#' 
#' # TBD
#' myPar <- tbdPar(c(4, 6))
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # GBCD
#' myPar <- gbcdPar(8, 2)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#' 
#' # ABCD
#' myPar <- abcdPar(8, 3)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#'
#' # BBCD
#' myPar <- bbcdPar(8, 5)
#' genSeq(myPar, 5)
#' genSeq(myPar)
#' 
#' # CHEN
#' myPar <- chenPar(8, 2, 0.667)
#' genSeq(myPar, 5)
#' genSeq(myPar)
#' 
#' @name generateRandomSequences
NULL


#' @rdname generateRandomSequences
#'
#' @export
setGeneric("genSeq", function(obj, r, seed) standardGeneric("genSeq"))


