##' Binding hyperSpec Objects
##' Two S3 functions \code{cbind.hyperSpec} and \code{rbind.hyperSpec} act as
##' an interfaces to \code{cbind2} and \code{rbind2} because neither
##' \code{\link[Matrix]{rBind}} and \code{\link[Matrix]{cBind}} nor S4 versions
##' of \code{cbind} and \code{rbind} do work at the moment.
##' 
##' While it is now possible to do S4 despatch on \code{\dots{}}, defining such
##' S4 methods for \code{cbind} and \code{rbind} breaks the binding of
##' \code{Matrix} objects. Therefore, two S3 methods \code{rbind.hyperSpec} and
##' \code{cbind.hyperSpec} are defined.
##' 
##' \code{bind} does the common work for both column- and row-wise binding.
##' 
##' @aliases bind
##' @param ... The \code{hyperSpec} objects to be combined.
##' 
##' Alternatively, \emph{one} list of \code{hyperSpec} objects can be given to
##'   \code{bind}.
##' @include paste.row.R
##' @param direction "r" or "c" to bind rows or columns
##' @return a \code{hyperSpec} object, possibly with different row order (for
##'   \code{bind ("c", \dots{})} and \code{cbind2}).
##' @note You might have to make sure that the objects either all have or all
##'   do not have rownames and/or colnames.
##' @author C. Beleites
##' @export 
##' @seealso \code{\link[Matrix]{rBind}}, \code{\link[Matrix]{cBind}}
##' \code{\link[methods]{rbind2}}, \code{\link[methods]{cbind2}}
##' \code{\link[base]{rbind}}, \code{\link[base]{cbind}}
##'
##' \code{\link{merge}} and \code{\link{collapse}} for combining objects that do not share spectra
##' or wavelengths, respectively.
##' @keywords methods manip
##' @examples
##' 
##' chondro
##' bind ("r", chondro, chondro)
##' rbind (chondro, chondro)
##' cbind (chondro, chondro)
##' bind ("r", list (chondro, chondro, chondro))
##' 
##' x <- chondro[,, 600 : 605]
##' x$a <- 1
##' x@@data <- x@@data[, sample (ncol (x), ncol (x))] # reorder columns
##' 
##' y <- chondro [nrow (chondro) : 1,, 1730 : 1750] # reorder rows
##' y$b <- 2
##' 
##' cbind2 (x, y) # works
##' 
##' y$y[3] <- 5
##' try (cbind2 (x, y)) # error
##' 
##' 
bind <- function (direction = stop ("direction ('c' or 'r') required"), ...){
  dots <- list (...)

  if ((length (dots) == 1) & is.list (dots [[1]]))
    dots <- dots[[1]]

  if (length (dots) == 0)
    NULL
  else if (length (dots) == 1){
    validObject (dots[[1]])
    dots[[1]]
  } else {                              # binding is actually needed.
    lapply (dots, chk.hy)
    lapply (dots, validObject)
    
    for (i in seq_along (dots) [-1]){
      dots[[1]] <- switch (direction,
                           c = cbind2 (dots[[1]], dots[[i]]),
                           r = rbind2 (dots[[1]], dots[[i]]),
                           stop ("direction must be either 'c' or 'r' for cbind",
                                 "and rbind, respectively.")
                           )
    }
    
    dots [[1]]
  }
}

##' \code{cbind2} binds the spectral matrices of two \code{hyperSpec} objects by column. All columns
##' besides \code{spc} with the same name in \code{x@@data} and \code{y@@data} must have the same
##' elements.  Rows are ordered before checking.
##' 
##' @aliases bind cbind.hyperSpec rbind.hyperSpec
##'   cbind2,hyperSpec,hyperSpec-method rbind2,hyperSpec,hyperSpec-method
##'   cbind2,hyperSpec,missing-method rbind2,hyperSpec,missing-method
##' @param x,y \code{hyperSpec} objects
##' @rdname bind
##' @export 
##' @aliases cbind.hyperSpec

##' 
cbind.hyperSpec <- function (...) bind ("c", ...)

##'
##' \code{rbind2} binds two \code{hyperSpec} objects by row. They need to have
##' the same columns.
##' 
##' @aliases  rbind.hyperSpec
##' @rdname bind
##' @export 
##' @aliases rbind.hyperSpec
rbind.hyperSpec <- function (...) bind ("r", ...)

##' @rdname bind
##' @export 
##' @aliases cbind2,hyperSpec,hyperSpec-method
setMethod ("cbind2", signature = signature (x = "hyperSpec", y = "hyperSpec"),
           function (x, y){
             validObject (x)
             validObject (y)

             cols <- match (colnames (x@data), colnames (y@data))
             cols <- colnames (y@data) [cols]
             cols <- cols [! is.na (cols)]
             cols <- cols [- match ("spc", cols)]

             if (length (cols) < 0){
               ord <- do.call (order, x@data[, cols, drop = FALSE])
               x@data <- x@data[ord, , drop = FALSE]

               ord <- do.call (order, y@data[, cols, drop = FALSE])
               y@data <- y@data[ord, , drop = FALSE]

               if (any (x@data[, cols, drop = FALSE] != y@data[, cols, drop = FALSE]))
                 stop ("hyperSpec objects must have the same data in columns",
                       "of the same name (except data$spc)")
             }

             ## for the spectra, multiple occurences of the same wavelength are O.K.
             x@data$spc <- cbind(x@data$spc, y@data$spc)
             .wl (x) <- c (x@wavelength, y@wavelength)

             ## cbind columns in y that are not in x
             cols <- is.na (match (colnames (y@data), colnames (x@data)))
             x@data <- cbind (x@data,
                              y@data[, cols, drop = FALSE])

             x
           }
           )

##' @rdname bind
##' @export 
##' @aliases cbind2,hyperSpec,missing-method
setMethod("cbind2", signature = signature (x = "hyperSpec", y = "missing"), function (x, y) x)

##' @rdname bind
##' @export 
##' @aliases  rbind2,hyperSpec,hyperSpec-method
setMethod("rbind2",
          signature = signature (x = "hyperSpec", y = "hyperSpec"),
          function (x, y) {
            validObject (x)
            validObject (y)

            if (! isTRUE (all.equal (x@wavelength, y@wavelength)))
              stop ("The wavelengths of the objects differ.\n",
                    "If they are not ordered, try 'orderwl'.")

            x@data <- rbind (x@data, y@data)

            x
          }
          )

##' @rdname bind
##' @export 
##' @aliases rbind2,hyperSpec,missing-method
setMethod ("rbind2", signature = signature (x = "hyperSpec", y = "missing"), function (x, y) x)

