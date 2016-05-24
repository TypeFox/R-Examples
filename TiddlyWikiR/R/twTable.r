##twTable.r
##2013-07-10 dmontaner@cipf.es
##TiddlyWikiR library

##' @name twTable-class
##' @docType class
##' 
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @aliases twTable
##' @aliases initialize,twTable-method
##' @aliases [,twTable-method
##' @aliases dat dat<-
##' @aliases ref,twTable-method  ref<-,twTable-method
##' @aliases color color<-
##' @aliases bgcolor bgcolor<-
##' @aliases includeRowNames includeRowNames<-
##' @aliases includeColNames includeColNames<-
##' @aliases rowref rowref<-
##' @aliases colref colref<-
##' @aliases align,twTable-method  align<-,twTable-method
##' @aliases title,twTable-method title<-
##' @aliases footer footer<-
##' @aliases sortable sortable<-
##' @aliases digits digits<-
##' @aliases dim,twTable-method
##' 
##' @keywords wiki table
##' @seealso \code{\link{twLink}} and \code{\link{twImage}}
##' 
##' @title A class to handle TiddlyWiki tables.
##'
##' @details Missing values are allowed for the entries of the table in the "dat" slot.
##' Missing data in in the "ref" or "color" slots are interpreted as no reference to be linked or not color to be used in the cell.
##'
##' If sortable = TRUE the plugin needs to be installed.
##' 
##' @section Usage:
##' new ("twTable", dat, ...)
##' 
##' twTable (dat, ref, color, ...)
##' 
##' @section Slots: 
##'  \describe{
##'    \item{\code{dat}:}{a data.frame containing the table values.}
##'    \item{\code{ref}:}{character matrix of the target or URL to be redirected to from each cell of the table.}
##'    \item{\code{color}:}{character matrix indicating the color of each cell of the table.}
##'    \item{\code{bgcolor}:}{character matrix indicating the background color of each cell of the table.}
##'    \item{\code{includeRowNames}:}{logical; if TRUE row names of "dat" are used as row names of the table.}
##'    \item{\code{includeColNames}:}{logical; if TRUE column names of "dat" are used as column names of the table.}
##'    \item{\code{rowref}:}{a character vectors of references to be linked from the row names of the table.}
##'    \item{\code{colref}:}{a character vectors of references to be linked from the column names of the table.}
##'    \item{\code{align}:}{a character vector indicating the alignment of each column of the table; may be missing or one of "l", "c" or "r".}
##'    \item{\code{title}:}{a title or header for the table.}
##'    \item{\code{footer}:}{a footer for the table.}
##'    \item{\code{sortable}:}{if TRUE the columns of the table may be sorted by clicking on their column names.}
##'    \item{\code{digits}:}{number of digits to be displayed in numerical, non integer columns. Default is 3.}
##'  }
##'
##' @examples
##' myData <- as.data.frame (matrix (rnorm (12), ncol = 4))
##' myTable <- new ("twTable", dat = myData)
##' wikify (myTable)
##' 
##' @export

setClass ("twTable",
          representation = representation (
            dat   = "data.frame",
            ref   = "matrix",
            color = "matrix",
            bgcolor = "matrix",
            ##
            ## rownames  = "character",
            ## colnames  = "character",
            ##
            includeRowNames = "logical",
            includeColNames = "logical",
            ##
            rowref = "character",
            colref = "character",
            ##
            align = "character",
            ##
            title  = "character",
            footer = "character",
            ##
            sortable = "logical",
            digits = "numeric"
            ))

##' @export
twTable <- function (...) {
  new ("twTable", ...)
}
################################################################################

##' @export
setMethod ("initialize", "twTable",
           function (.Object, 
                     dat, 
                     ref, 
                     color,
                     bgcolor,
                     ## rownames,
                     ## colnames, 
                     includeRowNames, 
                     includeColNames, 
                     rowref, 
                     colref, 
                     align, 
                     title, 
                     footer, 
                     sortable,
                     digits) { ##the variable has to be called .Object
             
             if (missing (dat)) {
               stop ("dat is missing with no default.")
             } else {
               if (is.matrix (dat)) {
                 dat <- as.data.frame (dat, stringsAsFactors = FALSE)
               }
               .Object@dat <- dat
             }

             ## ################################################################
             
             if (missing (ref)) {
               .Object@ref <- matrix (NA_character_ , nrow = nrow (dat), ncol = ncol (dat))
               dimnames (.Object@ref) <- dimnames (.Object@dat)
             } else {
               if (any (dim (dat) != dim (ref))) {
                 stop ('"dat" and "ref" have different dimensions.')
               } else {
                 .Object@ref <- ref
               }
             }

             if (missing (color)) {
               .Object@color <- matrix (NA_character_ , nrow = nrow (dat), ncol = ncol (dat))
               dimnames (.Object@color) <- dimnames (.Object@dat)
             } else {
               if (any (dim (dat) != dim (color))) {
                 stop ('"dat" and "color" have different dimensions.')
               } else {
                 .Object@color <- color
               }
             }

             if (missing (bgcolor)) {
               .Object@bgcolor <- matrix (NA_character_ , nrow = nrow (dat), ncol = ncol (dat))
               dimnames (.Object@bgcolor) <- dimnames (.Object@dat)
             } else {
               if (any (dim (dat) != dim (bgcolor))) {
                 stop ('"dat" and "bgcolor" have different dimensions.')
               } else {
                 .Object@bgcolor <- bgcolor
               }
             }

             ## ################################################################
             
             ## if (missing (rownames)) {
             ##   ##.Object@rownames <- rownames (dat) ##NOT WORKING
             ##   .Object@rownames <- dimnames (dat)[[1]]
             ## } else {
             ##   if (ncol (dat) != length (rownames)) {
             ##     stop ('"rownames" length does not match ncols in "dat".')
             ##   } else {
             ##     .Object@rownames <- rownames
             ##   }
             ## }
             
             ## if (missing (colnames)) {
             ##   ##.Object@colnames <- colnames (dat) ##NOT WORKING
             ##   .Object@colnames <- dimnames (dat)[[2]]
             ## } else {
             ##   if (ncol (dat) != length (colnames)) {
             ##     stop ('"colnames" length does not match ncols in "dat".')
             ##   } else {
             ##     .Object@colnames <- colnames
             ##   }
             ## }
             
             ## ################################################################

             if (missing (includeRowNames)) {
               .Object@includeRowNames <- TRUE
             } else {
               .Object@includeRowNames <- includeRowNames[1]
             }

             if (missing (includeColNames)) {
               .Object@includeColNames <- TRUE
             } else {
               .Object@includeColNames <- includeColNames[1]
             }

             ## ################################################################
             
             if (missing (rowref)) {
               .Object@rowref <- rep (NA_character_ , times = nrow (dat))
             } else {
               if (nrow (dat) != length (rowref)) {
                 stop ('"rowref" length does not match ncols in "dat".')
               } else {
                 .Object@rowref <- rowref
               }
             }
             
             if (missing (colref)) {
               .Object@colref <- rep (NA_character_ , times = ncol (dat))
             } else {
               if (ncol (dat) != length (colref)) {
                 stop ('"colref" length does not match ncols in "dat".')
               } else {
                 .Object@colref <- colref
               }
             }
             
             ## ################################################################

             ## Insert a space before cell content to right justify cell
             ## Insert a space after cell content to left justify cell
             ## Insert spaces before and after cell content to center justify cell
             ## Insert an exclamation mark (!) as the first non-space character of a cell to turn it into a header cell

             if (missing (align)) {
               .Object@align <- rep (NA_character_ , times = ncol (dat))
             } else {
               if (!all (align %in% c("l", "c", "r", NA))) {
                 stop ('"align" has to be one of: "l", "c", "r" or NA')
               }
               if (ncol (dat) != length (align)) {
                 stop ('"align" length does not match ncols in "dat".')
               } else {
                 .Object@align <- align
               }
             }

             ## ################################################################
             
             if (missing (title)) {
               .Object@title <- NA_character_
             } else {
               .Object@title <- title
             }
             
             if (missing (footer)) {
               .Object@footer <- NA_character_
             } else {
               .Object@footer <- footer
             }
             
             if (missing (sortable)) {
               .Object@sortable <- FALSE   ##make it depend on whether the plugging is installed or not
             } else {
               .Object@sortable <- sortable
             }

             if (missing (digits)) {
               .Object@digits <- 3
             } else {
               if (is.numeric (digits)) {
                 .Object@digits <- as.numeric (digits)
               } else {
                 stop ('"digits" has to be an integer.')
               }
             }
             
             ## ################################################################             
             
             return (.Object)
           })

################################################################################

## does not need @export
setValidity ("twTable", function (object) { ##the variable can be named other than "object"
  out <- NULL
  
  if (any (dim (object@dat) != dim (object@ref))) {
    out <- c (out, '"dat" and "ref" have different dimension.')
  }

  if (any (dim (object@dat) != dim (object@color))) {
    out <- c (out, '"dat" and "color" have different dimension.')
  }

  if (any (dim (object@dat) != dim (object@bgcolor))) {
    out <- c (out, '"dat" and "bgcolor" have different dimension.')
  }

  
  if (length (object@align) != ncol (object@dat)) {
    out <- c (out, '"align" length is not equal to the number of columns.')
  }
  
  if (!all (object@align %in% c("l", "c", "r", NA))) {
    out <- c (out, '"align" has to be one of: "l", "c", "r" or NA.')
  }

  if (!is.numeric (object@digits)) {
    out <- c (out, '"digits" has to be an integer.')
  }
  
  ##RETURN
  if (is.null (out)) {
    out <- TRUE
  }
  return (out)
})

################################################################################

##  Sub-setting operator
##' @export
setMethod ("[", "twTable", 
           function (x, i, j, ..., drop=TRUE) {
             
             if (!missing (i)) {
               if (is.character (i)) { ##use names
                 i <- match (i, rownames (x@dat))
                 i <- i[!is.na (i)]
                 print ("i texto")
               }
             }
             
             if (!missing (j)) {
               if (is.character (j)) { ##use names
                 j <- match (j, colnames (x@dat))
                 j <- j[!is.na (j)]
                 print ("j texto")
               }
             }

             ##logical index can be used already
             
             x@dat <- x@dat[i, j, drop = FALSE]
             x@ref <- x@ref[i, j, drop = FALSE]
             x@color <- x@color[i, j, drop = FALSE]
             x@bgcolor <- x@bgcolor[i, j, drop = FALSE]

             x@rowref <- x@rowref[i]
             x@colref <- x@colref[j]
             x@align  <- x@align[j]

             return (x)
           })

################################################################################
### ACCESSOR and REPLACEMENTS
################################################################################

##' @export
dat <- function (x) {
  x@dat
}

##' @export
'dat<-' <- function (x, value) {
  x@dat <- value
  return (x)
}

########################################

setMethod ("ref", "twTable", function (object) {
  return (object@ref)
})

setReplaceMethod ("ref", "twTable", function (object, value) {
  object@ref <- value
  return (object)
})

########################################

##' @export
color <- function (x) {
  x@color
}

##' @export
'color<-' <- function (x, value) {
  x@color <- value
  return (x)
}

########################################

##' @export
bgcolor <- function (x) {
  x@bgcolor
}

##' @export
'bgcolor<-' <- function (x, value) {
  x@bgcolor <- value
  return (x)
}

########################################

##' @export
includeRowNames <- function (x) {
  x@includeRowNames
}

##' @export
'includeRowNames<-' <- function (x, value) {
  x@includeRowNames <- value
  return (x)
}

########################################

##' @export
includeRowNames <- function (x) {
  x@includeRowNames
}

##' @export
'includeRowNames<-' <- function (x, value) {
  x@includeRowNames <- value
  return (x)
}

########################################

##' @export
rowref <- function (x) {
  x@rowref
}

##' @export
'rowref<-' <- function (x, value) {
  x@rowref <- value
  return (x)
}
########################################

##' @export
colref <- function (x) {
  x@colref
}

##' @export
'colref<-' <- function (x, value) {
  x@colref <- value
  return (x)
}

########################################

setMethod ("align", "twTable", function (object) {
  return (object@align)
})

setReplaceMethod ("align", "twTable", function (object, value) {
  object@align <- value
  return (object)
})

########################################

## title <- function (x) {  ##exists in the 'graphics' package
##   x@title
## }

##I am not sure this is the best way to define 'title' as generic and method function.
##Arguments: main, sub, xlab, ylab, line and outer
##need to be included as they are defined in the graphics::title which becomes the default method.
##Can I avoid graphics::title being the default so I do not need all these arguments in my method?
setMethod ("title", "twTable", function (object, main, sub, xlab, ylab, line, outer) {
  return (object@title)
})

##' @export
'title<-' <- function (x, value) {
  x@title <- value
  return (x)
}
########################################

##' @export
footer <- function (x) {
  x@footer
}

##' @export
'footer<-' <- function (x, value) {
  x@footer <- value
  return (x)
}
########################################

##' @export
sortable <- function (x) {
  x@sortable
}

##' @export
'sortable<-' <- function (x, value) {
  x@sortable <- value
  return (x)
}
########################################

##' @export
digits <- function (x) {
  x@digits
}

##' @export
'digits<-' <- function (x, value) {
  x@digits <- value
  return (x)
}

################################################################################

###Wikify

wikify.table <- function (object) {
  ## 1 link
  ## 2 rownames
  ## 3 align
  ## 4 color
  ## 5 background color
  ##SYNTAX:
  ##|color:#000099;TableTextHere|
  ##|bgcolor:#000099;TableTextHere|
  ##|bgcolor:#000099;color:#990000;TableTextHere|
  ##|bgcolor:#000099;color:#990000;[[TableTextHere|urlHere]]|
  ##the link overwrites the color
  
  ##FORMAT DATA
  dat     <- object@dat
  ref     <- object@ref
  color   <- object@color
  bgcolor <- object@bgcolor
  ##
  rowNames <- rownames (dat)
  colNames <- colnames (dat)
  Nrow <- nrow (dat)
  Ncol <- ncol (dat)

  ##rounding numerical
  es.numeric <- sapply (dat, class) == "numeric"
  if (any (es.numeric)) {
    dat[es.numeric] <- round (dat[es.numeric], digits = object@digits)
  }
  
  dat <- matrix (as.character (unlist (dat)), nrow = Nrow, ncol = Ncol) ##this tries to avoid the spaces when converting to text. FIND A NICER WAY TO DO IT.
  
  ##row and col names
  if (object@includeColNames) {
    dat     <- rbind (colnames (object@dat), dat)
    ref     <- rbind (object@colref,         ref)
    color   <- rbind (NA,                  color)
    bgcolor <- rbind (NA,                bgcolor)
    Nrow <- Nrow + 1
    if (object@includeRowNames) {
      dat     <- cbind (c ("", rowNames),        dat)
      ref     <- cbind (c (NA, object@rowref),   ref)
      color   <- cbind (NA,                    color)
      bgcolor <- cbind (NA,                  bgcolor)
      Ncol <- Ncol + 1
    }
  } else {
    if (object@includeRowNames) {
      dat     <- cbind (rowNames,        dat)
      ref     <- cbind (object@rowref,   ref)
      color   <- cbind (NA,            color)
      bgcolor <- cbind (NA,          bgcolor)
      Ncol <- Ncol + 1
    }
  }
  dimnames (dat) <- NULL

  
  ##include LINKS
  makelink <- function (x) {
    if (is.na(x[1])) {
      out <- ""
    } else {
      if (is.na(x[2])) {
        out <- x[1]
      } else {
        out <- paste ("[[", x[1], "|", x[2], "]]", sep = "")
      }
    }
  }

  lin <- apply (cbind (as.vector (dat), as.vector (ref)), 1, makelink)
  mat <- matrix (lin, nrow = Nrow, ncol = Ncol) ##this tries to avoid the spaces when converting to text. FIND A NICER WAY TO DO IT.

  ##format row names
  if (object@includeRowNames) {
    mat[,1] <- paste ("!", mat[,1], sep = "")
  }

  
  ##ALIGNMENT
  if (object@includeRowNames) {
    object@align <- c(NA, object@align)
  }
  align <- rep (object@align, each = Nrow)
  ##
  touse <- which (align %in% c("l", "c"))
  mat[touse] <- paste (mat[touse], "")
  ##
  touse <- which (align %in% c("r", "c"))
  mat[touse] <- paste ("", mat[touse])
  
  ##COLOR
  color <- as.character (color)
  con.color <- !is.na (color)
  mat[con.color] <- paste ("color:", color[con.color], ";", mat[con.color], sep = "")

  ##BGCOLOR
  bgcolor <- as.character (bgcolor)
  con.bgcolor <- !is.na (bgcolor)
  mat[con.bgcolor] <- paste ("bgcolor:", bgcolor[con.bgcolor], ";", mat[con.bgcolor], sep = "")

  
  ##WIKI FORMAT
  mat <- cbind ("", mat, "")
  mat <- apply (mat, 1, paste, collapse = "|")

  if (object@includeColNames) {
    mat[1] <- paste (mat[1], "h", sep = "")
  }

  
  ##SORTABLE
  if (object@sortable) {
    mat <- c ("|sortable|k", mat)
  }
  
  ##HEADER AND FOOTER
  if (!is.na (object@title)) {
    mat <- c (paste ("|", object@title, "|c", sep = ""), mat)
  }
  if (!is.na (object@footer)) {
    mat <- c (mat, paste ("|", object@footer, "|c", sep = ""))
  }

  ##newlines to separate tables
  mat <- c("", mat, "")
  
  ## OUTPUT
  return (mat)
}

###Wikify method
setMethod ("wikify", "twTable", wikify.table)


################################################################################

##' @export
setMethod ("dim", "twTable", function (x) dim (x@dat))
