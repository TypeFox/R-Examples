### xtable package
###
### Produce LaTeX and HTML tables from R objects.
###
### Copyright 2000-2013 David B. Dahl <dahl@stat.byu.edu>
###
### This file is part of the `xtable' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

"caption<-" <- function(x, value) UseMethod("caption<-")
"caption<-.xtable" <- function(x, value) {
  if (length(value) > 2)
    stop("\"caption\" must have length 1 or 2")
  attr(x, "caption") <- value
  return(x)
}

caption <- function(x, ...) UseMethod("caption")
caption.xtable <- function(x, ...) {
  return(attr(x, "caption", exact = TRUE))
}

"label<-" <- function(x, value) UseMethod("label<-")
"label<-.xtable" <- function(x, value) {
  if (length(value) > 1)
    stop("\"label\" must have length 1")
  attr(x, "label") <- value
  return(x)
}

label <- function(x, ...) UseMethod("label")
label.xtable <- function(x, ...) {
  return(attr(x, "label", exact = TRUE))
}

"align<-" <- function(x, value) UseMethod("align<-")

### Based on contribution from Jonathan Swinton <jonathan@swintons.net>
### in e-mail dated Wednesday, January 17, 2007
.alignStringToVector <- function(aString) {
  ## poor mans parsing - separating string of form "l{2in}llr|p{1in}c|{1in}"
  ## into "l{2in}" "l"  "l"  "r" "|" "p{1in}" "c" "|{1in}"
  aString.Align <- character(0);
  aString.Width <- character(0);
  wString <- aString
  while( nchar(wString) > 0) {
    aString.Align <- c(aString.Align, substr(wString, 1, 1))
    ## is it followed by a brace?
    thisWidth <- ""
    if ( nchar(wString) > 1 & substr(wString, 2, 2) == "{") {
      beforeNextBrace <- regexpr("[^\\]\\}", wString)
      if (beforeNextBrace <0 ) {
        stop("No closing } in align string")
      }
      thisWidth <- substr(wString, 2, beforeNextBrace + 1)
      wString <- substr(wString, beforeNextBrace + 2, nchar(wString))
    } else {
      wString <- substr(wString, 2, nchar(wString))
    }
    aString.Width <- c(aString.Width, thisWidth)
  }

  alignAllowed <- c("l","r","p","c","|","X")

  if (any( !(aString.Align %in% alignAllowed))) {
    warning("Nonstandard alignments in align string")
  }
  res <- paste(aString.Align, aString.Width, sep = "")
  res
}
###.alignStringToVector ("l{2in}llr|p{1in}c|{1in}")
###.alignStringToVector ("l{2in}llr|p{1in}c|")
### latex syntax error, but gives wrong alignment
###.alignStringToVector ("{2in}llr|p{1in}c|")
###.alignStringToVector("llllp{3cm}")

"align<-.xtable" <- function(x, value) {
### Based on contribution from Benno <puetz@mpipsykl.mpg.de>
### in e-mail dated Wednesday, December 01, 2004
### Based on contribution from Jonathan Swinton <jonathan@swintons.net>
### in e-mail dated Wednesday, January 17, 2007
  ## cat("%", value, "\n")
  if ( (!is.null(value)) && ( is.character(value) ) &&
       ( length(value) == 1 ) && ( nchar(value) > 1 ) ) {
        value <- .alignStringToVector(value)
  }
  ## That should have checked we had only lrcp|
  ## but what if the "if statement" is false?
  ## For simplicity, deleting check present in version 1.4-2 and earlier.
  c.value <- if (any(!is.na(match(value, "|")))) {
               value[-which(value == '|')]
             } else {
               value
             }
  if (length(c.value) != ncol(x) + 1)
    stop(paste("\"align\" must have length equal to",
               ncol(x) + 1, "( ncol(x) + 1 )"))
  attr(x, "align") <- value
  return(x)
}

align <- function(x, ...) UseMethod("align")
align.xtable <- function(x, ...) {
  return(attr(x, "align", exact = TRUE))
}

"digits<-" <- function(x, value) UseMethod("digits<-")
"digits<-.xtable" <- function(x, value) {
  if( is.matrix( value ) ) {
    if( ncol( value ) != ncol(x) + 1 || nrow( value ) != nrow(x) ) {
      stop( "if argument 'digits' is a matrix, it must have columns equal",
           " to ", ncol(x) + 1, " ( ncol(x) + 1 ) and rows equal to ", nrow(x),
           " ( nrow( x )" )
    }
  } else {
    if( length(value) == 1 ) { value <- rep(value, ncol(x) + 1) }
    if( length( value ) > 1 & length( value ) != ncol(x) + 1 ) {
      stop( "if argument 'digits' is a vector of length more than one, it must have length equal",
           " to ", ncol(x) + 1, " ( ncol(x) + 1 )" )
    }
  }
  if (!is.numeric(value))
    stop("\"digits\" must be numeric")
  attr(x, "digits") <- value
  return(x)
}

digits <- function(x, ...) UseMethod("digits")
digits.xtable <- function(x, ...) {
  return(attr(x, "digits", exact = TRUE))
}

"display<-" <- function(x, value) UseMethod("display<-")
"display<-.xtable" <- function(x, value) {
  if (length(value) != ncol(x) + 1)
    stop(paste("\"display\" must have length equal to",
               ncol(x) + 1, "( ncol(x) + 1 )"))
  if (!all(!is.na(match(value, c("d","f","e","E","g","G","fg","s")))))
    stop("\"display\" must be in {\"d\",\"f\",\"e\",\"E\",\"g\",\"G\", \"fg\", \"s\"}")
  attr(x, "display") <- value
  return(x)
}

display <- function(x, ...) UseMethod("display")
display.xtable <- function(x, ...) {
  return(attr(x, "display", exact = TRUE))
}

