###############################################
# --------------------------------------------#
# Class randSeq                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

validateRandSeq <- function(object) {
  errors <- character()
  groups <- object@groups
  K <- object@K


  if (!(length(groups) == K)) {
    msg <- paste("Length of groups is ", length(groups), ". Should have length ", K,
                 ".", sep = "")
    errors <- c(errors, msg)
  }
  
  if (.hasSlot(object, "seed")) {
    seed <- object@seed
    if (length(seed) != 1) {  
      warning(paste("Length of seed is ", length(seed), ". First argument ", seed[1] ,
                    " is used.", sep = ""))

    }
    if (!(round(seed[1]) == seed[1])) {  
      warning(paste("First argument of seed is ", seed[1], ". Used seed was", 
                    integer(seed[1]), ".", sep = ""))
      
    }
  }
  
  if (length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for randSeq
# --------------------------------------------

#' @title An S4 Class for the representation of  randomization sequences
#' 
#' @description This set of classes provides functionality of storing randomization
#' sequences of different randomization procedures along with the parameters 
#' representing the design.
#' 
#' @slot N total number of patients included in the trial
#' @slot M matrix containing randomization sequences of length \code{N} in its
#' rows.
#' @slot K number of treatment groups
#' @slot groups character string of length K defining the names of the treatment groups
setClass("randSeq",
         slots = c(M = "matrix", N = "numeric", K = "numeric",
                   ratio = "numeric", groups = "character"),
          validity = validateRandSeq)

# --------------------------------------------
# Class definition for rRandSeq
# --------------------------------------------

# @title An S4 Class for the representation of  randomization sequences generated at random.
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of different randomization procedures along with the parameters 
# representing the design 
# 
# @slot seed integer specifying the seed for the generation of randomization sequences
setClass("rRandSeq",
         slots = c(seed = "numeric"),
		 contains = "randSeq",
         validity = validateRandSeq)

# --------------------------------------------
# Accesssor functions for randSeq
# --------------------------------------------

#' Method defining the $ operator for the randSeq class
#' 
#' @inheritParams overview
setMethod("$", "randSeq",
          function(x, name) slot(x, name))

		  
#' Function returning the allocation seed slot of an object
#'
#' Returns the seed that was either generated at random or user specified.
#' The seed can be specified for any random operation e.g. genSeq.
#'
#' @inheritParams overview
seed <- function(obj) {
  if (.hasSlot(obj, "seed")) obj@seed
  else stop("Object has no slot named seed.") 
}

#' Accessor function for the randomization list 
#'
#' Get the randomization list coded in its groups.
#'
#' @inheritParams overview 
#'
#' @examples 
#' myPar <- bsdPar(10, 2)
#' M <- genSeq(myPar, 2)
#' getRandList(M)
#'
#' @name getRandomizationList
#'
#' @export
getRandList <- function(obj) {
  if (.hasSlot(obj, "M")) {
	sequences1 <- sequences2 <- obj@M        
    for(i in 1:obj@K) {
      sequences1[sequences2 == i-1] <- obj@groups[i]
    }
    sequences1
  }	
  else stop("Object has no slot named M.") 
}


# --------------------------------------------
# Show function for randSeq
# --------------------------------------------

setMethod("show", "randSeq", function(object) {
  # headline
  cat("\nObject of class \"", class(object)[1],"\"\n\n", sep = "")
  # crop the method from the class name of the randPar object
  cat("design =", getDesign(object), "\n") 
  # iterate through all slots of the object
  names <- slotNames(object)
  names <- names[!(names == "M")] 
  if (K(object) == 2) names <- names[!(names %in% "K")]
  if (all(ratio(object) == rep(1, K(object)))) {
    names <- names[!(names %in% "ratio")]
  }
  for(name in names) {
    cat(name, "=", slot(object, name), "\n")
  }  
  # The matrix M is printed seperately dependent on its size.
  print.matrix <- function(m) {
    write.table(format(m, justify = "left"),
            row.names = T, col.names = F, quote = F)
  }

  if (nrow(object@M) %in% 2:3) {
    sequences <- getRandList(object)
    cat("\nThe sequences M: \n\n")
    if (ncol(sequences) < 11) {
      print(sequences)
    } else {
      print(cbind(sequences[ , 1:10], "..."))
    }
  } else if (nrow(object@M) == 1) {
     sequences <- getRandList(object)
     cat("\nThe sequence M: \n\n")
     if (ncol(sequences) < 11) {
       print(sequences)
     } else {
       sequences <- t(matrix(sequences[1,1:11]))
       sequences[1,11] <- "..."
       print(sequences)
     }
  } else {
      cat("\nThe first 3 of", nrow(object@M), "sequences of M: \n\n")
      object@M <- object@M[1:4, ]
      sequences <- getRandList(object)
      if (ncol(sequences) < 11) {
        print(sequences[1:3, ])
        cat("...")
      } else {
        print(cbind(sequences[1:3, 1:10], "..."))
        cat("...")
      }
  }
  
  cat("\n") 
})


# --------------------------------------------
# Generic functions for randSeq
# --------------------------------------------

#' Theoretical probability for randomization sequences
#'
#' Calculate theoretical probability for observed randomization sequences
#'
#' @aliases getProbabilities calculateProbabilities calcProb
#' 
#' @param obj object of a class inheriting from randSeq. Formal representation 
#' of a randomization sequences together with the parameters that belong to
#' the procedure that generated the sequences.
#'
#' @examples 
#' myPar <- bsdPar(10, 2)
#' M <- genSeq(myPar, 2)
#' getProb(M)
#' 
#' # All Sequences
#' par <- pbrPar(bc=c(2,2))
#' refSet <- getAllSeq(myPar)
#' probs <- getProb(refSet)
#' 
#' # Sequences with probabilities
#' cbind(probs, refSet$M)
#' 
#' @name getProbabilities
NULL

#' @rdname getProbabilities
#'
#' @export
setGeneric("getProb", function(obj) standardGeneric("getProb"))


#' Sequence plotting
#'
#' Plot all randomization sequences of a randSeq object
#' 
#' @param sequences object of type randSeq
#' @param emph integer indicating which sequence should be highlighted in blue.
#' @param plotAllSeq logical. If \code{plotAllSeq=TRUE}, the complete set of 
#' randomization sequences will be plotted in light gray.
#' @param rs vector of a randomization sequence that should be highlighted.
#' 
#' @export
plotSeq <- function(sequences, plotAllSeq = FALSE, emph = NA, rs = NA){ 
  N <- N(sequences)
  
  plot.new()
  plot.window(xlim = c(0, N), ylim = c(-N, N))
  abline(a = 0, b = 0, col="lightgray")
  axis(1)
  axis(2)
  axis(4)
  #title(main = "Randomization Sequences")
  title(ylab = "Difference in group size")
  title(xlab = expression(paste("Patient ", i)))
  box()
  
  if (plotAllSeq) {
    for (i in 0:(N-1)) {
      for (j in seq(-i,i,2)) {
        lines(c(i, i+1), c(j, j+1), type = "b", col = "lightgray")
        lines(c(i, i+1), c(j, j-1), type = "b", col = "lightgray")
      }
    }
  }
  
  numberOfSequences <- nrow(sequences@M)
  if (!is.na(emph)) stopifnot(emph < numberOfSequences, 0 < emph)
  
  for (i in 1:numberOfSequences) {
    lines(0:N, cumsum(c(0,2*(sequences@M)[i,]-1)), type = "b")
  }
  
  if(!is.na(emph)){
    lines(0:N, cumsum(c(0,2*sequences@M[emph,]-1)), type = "b",lwd = 2, col = "cornflowerblue")
  }
  else if (!anyNA(rs)) {
    f<-c("black", "red")
    D<-c(0, cumsum(2*as.numeric(rs)-1)) # contains the random walk
    lines(0:N, D, type = "b", lwd=2)
  }
}


