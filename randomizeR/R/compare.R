#' @include issue.R
#' @include randSeq.R
#' @include chronBias.R
#' @include corGuess.R
#' @include util.R
#' @include endpoint.R
#' @include assess.R
NULL

###############################################
# --------------------------------------------#
# Class comparison                            #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the comparison class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateComparison <- function(object) {
  errors <- character()
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for comparison
# --------------------------------------------

# Randomization paramters generic
setClass("comparison",
         slots = c(S = "data.frame", L = "list"),
         validity = validateComparison)


# --------------------------------------------
# Accesssor functions for comparison
# --------------------------------------------

#' Method defining the $ operator for the assessemnt class
#' 
#' @inheritParams overview
setMethod("$", "comparison",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for comparison
# --------------------------------------------

setMethod("show", "comparison", function(object) {
  validObject(object)
  # headline
  cat("\nComparison for ", colnames(object$L[[1]]$D)[3],"\n\n", sep = "")
  print(round(object@S, digits = 3))
  cat("\n") 
})


# --------------------------------------------
# Generic functions for using objects of type issue and randSeq
# --------------------------------------------

#' Comparison of randomization procedures
#'
#' Compares randomization procedures based on a specified issue 
#' in clinical trials.
#'
#' @param issue object of class \code{issue}.
#' @param endp object of class \code{endpoint}, or \code{missing}.
#' @param ... at least one object of class \code{randSeq} or a list of 
#' objects of class \code{randSeq}.
#'
#' @details
#' Randomization procedures behave differently with respect to issues
#' like selection bias, chronological bias, or loss in power estimation.
#' The \code{compare} function evaluates the behaviour of randomization 
#' procedures with respect to one issue. 
#' Its first argument should represent one of the implemented 
#' \code{\link{issues}}.
#' The second argument should be any number of objects of the class
#' \code{randSeq}. These objects represent the randomization procedures
#' for the planned comparison. 
#' The last argument \code{endp} may be provided if 
#' the assessment should take the distribution of the treamtent groups
#' into account, e.g. for power evaluation.
#'
#' @examples 
#' # compare Random Allocation Rule and Big Stick for N = 4
#' # with respect to the correct guesses
#' RAR <- getAllSeq(rarPar(4))
#' BSD <- getAllSeq(bsdPar(4, mti = 2))
#' corGuess <- corGuess("CS")
#' (comp <- compare(corGuess, RAR, BSD))
#' plot(comp)
#'
#' # compare the same procedures with respect to selection bias
#' endp <- normEndp(c(2, 2), c(1, 1))
#' selBias <- selBias("CS", 4, "exact")
#' (comp <- compare(selBias, RAR, BSD, endp = endp))
#' plot(comp)
#'
#' @return
#' \code{S4} object of class \code{comparison} summarizing the comparison of the 
#' randomization procedures.
#'
#' @seealso Representation of randomization procedures: \code{\link{randPar}}
#' @seealso Generation of randomization sequences: \code{\link{genSeq}}
#' @seealso \code{\link{issues}} for the assessment of randomization sequences
#' 
#' @name compare
NULL

#' @rdname compare
#'
#' @export
setGeneric("compare", function(issue, ..., endp) standardGeneric("compare"))


#' Generic plotting of comparison objects
#' 
#' @param x object of class \code{comparison}.
#' @param y character \code{"boxplot"}, or \code{"violin"}, or \code{"missing"}.
#' @param ... \code{"missing"}
#'
#' @details
#' Creates a box- or violinplot of an object  \code{x} of the class \code{comparison}. 
#' 
#' @return 
#' A plot created with the additional package \code{ggplot2}.
#' 
#' @examples 
#' # compare Random Allocation Rule and Big Stick for N = 4
#' # with respect to the correct guesses
#' RAR <- getAllSeq(rarPar(4))
#' BSD <- getAllSeq(bsdPar(4, mti = 2))
#' corGuess <- corGuess("CS")
#' comp <- compare(corGuess, RAR, BSD)
#' plot(comp)
#' 
#' @seealso \code{\link{compare}} for creating \code{S4} objects of the class \code{comparison}
#' 
#' @name plot
NULL

#' @rdname plot
#'
#' @export
setGeneric("plot")

# --------------------------------------------
# Methods for assessment
# --------------------------------------------

#' @rdname compare
setMethod("compare", signature(issue = "issue", endp = "missing"),
          function(issue, ...) {
            R <- list(...)
            if (length(R) == 1 && is.list(R[[1]])) {
              R <- c(...)
            }
            
            if (!all(sapply(R, function(x)  is(x, "randSeq"))))
              stop("Not all ... objects of class randSeq.")
            
            if (!all(sapply(R, function(x)  identical(x@K, 2))))
              stop("Not all ... objects have K = 2.")
            
            if (!all(sapply(R, function(x)  identical(x@ratio, c(1, 1)))))
              stop("Not all ... objects have ratio = c(1, 1).")
            
            output <- lapply(R, function(r) {
             assess(r, issue)  
            })
            
            names <- lapply(output, function(x) {
              x$design
            })
            
            S <- do.call("cbind", lapply(output, function(x) {
              summary(x)
            })) 
            colnames(S) <- names
            
            new("comparison", S = data.frame(S), L = output)  
          }
)

#' @rdname compare
setMethod("compare", signature(issue = "issue", endp = "endpoint"),
          function(issue, ..., endp) {
            R <- list(...)
            if (length(R) == 1 && is.list(R[[1]])) {
              R <- c(...)
            }
            
            if (!all(sapply(R, function(x)  is(x, "randSeq"))))
              stop("Not all ... objects of class randSeq.")
            
            if (!all(sapply(R, function(x)  identical(x@K, 2))))
              stop("Not all ... objects have K = 2.")
            
            if (!all(sapply(R, function(x)  identical(x@ratio, c(1, 1)))))
              stop("Not all ... objects have ratio = c(1, 1).")
            
            output <- lapply(R, function(r) {
              assess(r, issue, endp = endp)  
            })
            
            names <- lapply(output, function(x) {
              x$design
            })
            
            S <- do.call("cbind", lapply(output, function(x) {
              summary(x)
            })) 
            colnames(S) <- names
            
            new("comparison", S = data.frame(S), L = output)  
          }
)


# --------------------------------------------
# Plot function for comparison
# --------------------------------------------

#' @rdname plot
setMethod("plot", c("comparison","character"), function(x,y) {
  Probability <- design <- NULL
  obj <- do.call("rbind", lapply(x@L, function(z) {
    D <- z@D
    if (sum(!duplicated(D[,3])) == 1 && y == "violin") {
       stop(paste("Violin plot not possible. All entries for", colnames(D)[3], "in", z@design, "are the same."))
    }

    D$Sequence <- NULL
    colnames(D)[1] <- "Probability"
    D$design <- z@design
    D
  }
  ))
  xName <- colnames(obj)[2]
  colnames(obj)[2]<- "x"
  if (is.logical(obj$x)) {
    warning("Issue is measured on a logical scale. Method should be adjusted to exact instead of sim.")
  }
    
  g <- ggplot(obj, aes(x = design, y = x, weight = Probability))
  if (y == "boxplot") {
    g <- g + geom_boxplot(aes(fill = design))
  }
  else { 
    g <- g + geom_violin(aes(fill = design))
  }
  g <- g + theme(legend.position = "none")
  g <- g + scale_y_continuous(name = xName)
  g <- g + scale_x_discrete(name = "Randomization Procedures")
  plot(g)
})

#' @rdname plot
setMethod("plot", c("comparison","missing"), function(x,y) {
  Probability <- design <- NULL
  obj <- do.call("rbind", lapply(x@L, function(z) {
    D <- z@D
    if (sum(!duplicated(D[,3])) == 1 && y == "violin") {
      stop(paste("Violin plot not possible. All entries for", colnames(D)[3], "in", z@design, "are the same."))
    }
    
    D$Sequence <- NULL
    colnames(D)[1] <- "Probability"
    D$design <- z@design
    D
  }
  ))
  xName <- colnames(obj)[2]
  colnames(obj)[2]<- "x"
  if (is.logical(obj$x)) {
    warning("Issue is measured on a logical scale. Method should be adjusted to exact instead of sim.")
  }
  
  g <- ggplot(obj, aes(x = design, y = x, weight = Probability))
  g <- g + geom_violin(aes(fill = design))
  g <- g + theme(legend.position = "none")
  g <- g + scale_y_continuous(name = xName)
  g <- g + scale_x_discrete(name = "Randomization Procedures")
  plot(g)
})
