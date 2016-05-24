#' Tools for functional programming in R
#'
#' This package contains a collection of functions that facilitate modeling 
#' of data using a functional programming paradigm. The idea is that 
#' using tools that are more closely connected with the 
#' idioms of mathematics will make it easier to map 
#' the mathematical model to the software model. 
#'
#' \tabular{ll}{
#' Package: \tab lambda.tools\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.7\cr
#' Date: \tab 2015-06-09\cr
#' License: \tab LGPL-3\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @section Details:
#' Functional programming concepts start with functions as the foundation.
#' Higher-order functions provide generalized machinery for operating
#' on data in an element-wise manner. Lambda.tools includes idiomatic
#' versions of the canonical higher-order functions, such as map and fold
#' for data structures common in R.
#' In most languages the semantics are
#' limited to single-element iterations. In R it is common to work with
#' panel data or sliding windows, so \code{lambda.tools} introduces block
#' and range semantics to support these concepts, respectively. 
#' Hence \code{lambda.tools} defines \code{mapblock} and \code{maprange}
#' and similar functions for \code{fold}. 
#'
#' \subsection{Block operations}{
#' The semantics of a block operation is that regular, continguous chunks of
#' data are passed to the function. Suppose a vector \code{x} has 12 elements.
#' Performing a mapblock operation with window of length 3 applies the 
#' specified function to the following sub-vectors: \code{x[1:3]},
#' \code{x[4:6]}, \code{x[7:9]}, \code{x[10:12]}.
#' This is useful for processing any vector or 
#' list produced by a function that returns a regular length output.
#' 
#' Note that if the original sequence is not an integer multiple of the
#' window length, the last sub-vector will not have the same length as
#' the preceding sub-vectors.
#' }
#' 
#' \subsection{Range operations}{
#' While block operations use adjacent sub-vectors, range operations
#' use overlapping sub-vectors. This process is analogous to a 
#' sliding window, where the index increments by one as opposed to
#' by the window size. For the same vector \code{x}, a maprange operation
#' with window of length 3 produces the following sub-vectors as
#' arguments: \code{x[1:3]}, \code{x[2:4]}, \code{x[3:6]}, ..., \code{x[10:12]}.
#' 
#' An example of a range operation is generating n-grams from a text
#' document. Suppose a vector \code{v} contains a sequence of words. Then
#' \code{maprange(v, 2, function(x) paste(x, collapse=' '))} creates bigrams.
#' }
#' 
#' \subsection{Two-dimensional operations}{
#' Typically map and fold operate on 1-dimensional data structures,
#' but in R operations can also be applied on 2-dimensional data structures.
#' For example, the \code{apply} function works in this manner where the
#' \code{MARGIN} argument defines whether iteration operates on rows versus
#' columns. 
#' Hence \code{lambda.tools} introduces 2-dimensional versions of
#' these functions.  For simplicity, the 2-dimensional variants of 
#' map and fold only operate along columns.
#' To operate along rows requires transposing the data structure.
#' 
#' Consider the following code that applies multiple
#' rotations to a collection of points.
#' 
#'   ps <- t(matrix(c(0,0, 4,0, 2,4), nrow=2))
#'   rt <- matrix(c(cos(pi),-sin(pi),sin(pi),cos(pi), 
#'     cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)), nrow=2)
#'   mapblock(rt, 2, function(x) ps %*% x)
#' 
#' The result is a 6x2 matrix that is the union of the two rotation
#' operations.
#' }
#'
#' \subsection{Other goodies}{
#' Other functions included are functions to manipulate sequences,
#' such as \code{pad} a sequence to a specified length, \code{chomp}
#' the head and tail off a vector, \code{slice} a sequence into
#' two pieces based on an expression. The \code{partition} function
#' is similar, while \code{quantize} and \code{confine} transform data
#' to fit specific ranges.
#'
#' Logical functions such as \code{onlyif} and \code{use_default}
#' eliminate the need for conditional blocks, which
#' can streamline code and remove the risk of poorly scoped variables.
#' }
#'
#' @name lambda.tools-package
#' @aliases lambda.tools-package lambda.tools
#' @docType package
#' @exportPattern "^[^\\.]"
#' @import lambda.r
#' @author Brian Lee Yung Rowe <r@@zatonovo.com>
#' @references Rowe, Brian Lee Yung.
#' Modeling Data With Functional Programming In R. Chapman & Hall/CRC Press.
#' Forthcoming.
#' @seealso \code{\link{map}} \code{\link{fold}} \code{\link{samplerange}}
#'  \code{\link{slice}} \code{\link{onlyif}}
#'  \code{\link{quantize}} \code{\link{partition}}
#'  \code{\link{lambda.r}}
#' @keywords package attribute logic
NULL
