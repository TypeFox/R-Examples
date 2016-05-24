##' Conversion of a hyperSpec object into a data.frame or matrix
##' \code{as.data.frame} returns \code{x@@data} (as data.frame) \code{as.matrix}
##' returns the spectra matrix \code{x@@data$spc} as matrix
##'
##' @rdname asdataframe
##' @name as.data.frame
##' @aliases as.data.frame as.data.frame,hyperSpec-method as.data.frame.hyperSpec
##' @docType methods
##' @param x a \code{hyperSpec} object
##' @param row.names if \code{TRUE}, a column \code{.row} is created containing row names or row
##' indices if no rownames are set. If character vector, the rownames are set accordingly.
##' @param optional ignored
##' @return \code{x@@data} and \code{x@@data$spc} (== \code{x$spc} == \code{x [[]]}), respectively.
##' @author C. Beleites
##' @method as.data.frame hyperSpec
##' @export
##' @seealso \code{\link[base]{as.data.frame}} 
##' @keywords methods
##' @examples
##' 
##' as.data.frame (chondro [1:3,, 600 ~ 620])
##' as.matrix (chondro [1:3,, 600 ~ 620])
##' lm (c ~ spc, data = flu [,,450])

as.data.frame.hyperSpec <- function (x, row.names = TRUE, optional =  NULL, ...){
  validObject (x)

  x <- x@data
  
  if (!is.null (row.names))
    if (isTRUE (row.names)){
      if (is.null (rownames (x)))
        x$.row <- seq_len (nrow (x))
      else
        x$.row <- rownames (x)
    } else
      rownames (x) <- row.names
  
  x
}
##' @method as.matrix hyperSpec
##' @rdname asdataframe
##' @param ... ignored 
##' @aliases as.matrix as.matrix,hyperSpec-method
##' @export
##' @seealso and \code{\link[base]{as.matrix}}
##' 
##' \code{\link[hyperSpec:extractreplace]{[}} for a shortcut to \code{as.matrix}
as.matrix.hyperSpec <- function (x, ...){
  validObject (x)

  unclass (x@data$spc)                  # remove the AsIs
}


##' \code{as.wide.df} converts the spectra matrix to a data.frame. The extra
##' data together with this data is returned. The column names of the spectra
##' matrix are retained (if they are numbers, without preceeding letters).
##' 
##' @rdname asdataframe
##' @aliases  as.wide.df
##' @export
##' @return 
##' 
##' \code{as.wide.df} returns a data.frame that consists of the extra data and
##'   the spectra matrix converted to a data.frame. The spectra matrix is
##'   expanded \emph{in place}.
##' @examples
##' 
##' as.wide.df (chondro [1:5,, 600 ~ 610])
##' summary (as.wide.df (chondro [1:5,, 600 ~ 610]))

as.wide.df <- function (x) {
  chk.hy (x)
  validObject (x)

  ispc <- match ("spc", colnames (x@data))

  ## logical indexing creates n by 0 data.frame that can be cbound, thus
  ## avoiding trouble with empty or 0 x 0 data.frames:
  before <- seq_len (ncol (x@data)) < ispc
  after <- seq_len (ncol (x@data)) > ispc

  ## colnames should be preserved

  cols <- c (colnames (x@data)  [before],
             colnames (x@data$spc),
             colnames (x@data) [after])

  x <- cbind (x@data [, before],
              as.data.frame (unclass (x@data [, ispc])),
              x@data [, after])
  colnames (x) <- cols
  x
}


##' \code{as.long.df} returns a long-format data.frame.
##' 
##' The data.frame returned by \code{as.long.df} is guaranteed to have columns
##' \code{spc} and \code{.wavelength}. If \code{nwl (x) == 0} these columns
##' will be \code{NA}.
##' 
##' @rdname asdataframe
##' @aliases as.long.df 
##' @param rownames should the rownames be in column \code{.rownames} of the
##'   long-format data.frame?
##' @param wl.factor should the wavelengths be returned as a factor (instead of
##'   numeric)?
##' @param na.rm if \code{TRUE}, rows where spc is not \code{NA} are deleted.
##' @export
##' @return \code{as.long.df} returns the stacked or molten version of \code{x@@data}. The
##'   wavelengths are in column \code{.wavelength}.
##' @seealso
##'  
##' \code{\link[utils]{stack}} and \code{\link[reshape]{melt}} or \code{\link[reshape2]{melt}} for
##' other functions producing long-format data.frames.
##' @examples
##' 
##' as.long.df (flu [,, 405 ~ 410])
##' summary (as.long.df (flu [,, 405 ~ 410]))
##' summary (as.long.df (flu [,, 405 ~ 410], rownames = TRUE))
##' summary (as.long.df (flu [,, 405 ~ 410], wl.factor = TRUE))
##' 
as.long.df <- function (x, rownames = FALSE, wl.factor = FALSE, na.rm = TRUE) {
  chk.hy (x)
  validObject (x)

  ispc <- match ("spc", colnames (x@data))

  if (nwl (x) == 0) {
    tmp <- cbind (data.frame (.wavelength = rep (NA, nrow (x)),
                              spc = rep (NA, nrow (x))),
                  x@data [, -ispc, drop = FALSE])
  } else {
    tmp <- x@data [rep (row.seq (x), nwl (x)), -ispc, drop = FALSE]

    tmp <- cbind (data.frame (.wavelength = rep (x@wavelength, each = nrow (x)),
                              spc = as.numeric (x [[]])),
                  tmp)
    if (wl.factor){
      tmp$.wavelength <- as.factor (tmp$.wavelength)
      wl <- colnames (x@data$spc)       # there may be a fancily formatted version in the column
                                        # names
      if (is.null (wl))
        wl <- x@wavelength              # if not, use the wavelength vector
      
      levels (tmp$.wavelength) <- wl
    }
  }

  if (rownames)
    tmp <- data.frame (.rownames = as.factor (rep (rownames (x),
                         length.out = nrow (tmp))),
                       tmp)

  if (na.rm)
    tmp <- tmp [!is.na (tmp$spc), ]
  
  tmp
}

##' 
##' \code{as.t.df} produces a 'transposed' data.frame with columns containing the spectra.
##' @rdname asdataframe
##' @aliases as.t.df
##' @return \code{as.t.df} returns a data.frame similar to \code{as.long.df}, but each
##'   spectrum in its own column. This is useful for exporting summary spectra,
##'   see the example.
##' @export
##' @examples
##' df <- as.t.df (apply (chondro, 2, mean_pm_sd))
##' head (df)
##' 
##' if (require (ggplot2)){
##'   ggplot (df, aes (x = .wavelength)) +
##'     geom_ribbon (aes (ymin = mean.minus.sd, ymax = mean.plus.sd),
##'       fill = "#00000040") +
##'     geom_line (aes (y = mean))
##' }
as.t.df <- function (x) {
  chk.hy (x)
  validObject (x)

  df <- as.data.frame (t (unclass (x@data$spc)))
  colnames (df) <- rownames (x@data)
  
  cbind (.wavelength = x@wavelength, df)
}
