#' Binary Logic GNU R Package

#' @description 
#' This package contains the \strong{\link{binary}} S3 class.
#' A data object can be instantiated to store a binary number(Base2).
#'
#' It can be used to convert, negate, shift or rotate the binary number.
#' (switchEndianess, bytesNeeded, binaryPrefix, fillUpToByte).
#' 
#' Binary operators:
#' \itemize{
#'  \item \strong{==} , \strong{!=} , \strong{<} , \strong{<=} , \strong{>} , \strong{>=}
#'  \item \strong{+} , \strong{-} , \strong{^}, \strong{*}
#'  \item \strong{&} , \strong{|} , \strong{xor} (Logical Operator. Bitwise operation. The smaller vector is added up with zeros)
#'  \item \strong{!} (Indicates logical negation (NOT). Bitwise Operations)
#' }
#' 
#' binaryLogic functions:
#' \itemize{
#'  \item \link{shiftLeft}(binary) , \link{shiftRight}(binary)
#'  \item \link{rotate}(binary)
#'  \item \link{negate}(binary)
#'  \item \link{switchEndianess}(binary)
#' }
#'
#' Additional function:
#' \itemize{
#'  \item \link{fillUpToByte}, \link{fillUpToBit}
#'  \item \link{bytesNeeded}
#'  \item \link{binaryPrefix}
#'  \item  \link{byte}
#' }
#'
#' @details
#' This \strong{\link{binary}} class is just not that great at heavy number crunching, but it brings some benefits.
#' Especially if you like to work using vectors in R. It inherits from the \emph{logical} class.
#' Some function from package \pkg{binaryLogic} can be applied to \emph{logical} vectors. Such as shift or rotate (see help).
#' 
#' The internal structure looks like this:
#' 
#' \code{structure(c(TRUE, FALSE), class = c("binary", "logical"), signed = FALSE, littleEndian = FALSE)}
#' 
#' It is composed of a \emph{logical} vector and several attributes. This structure shows a big endian number, 
#' it corresponds to the value = 2 (Base10).
#' @docType package
#' @name binaryLogic
NULL