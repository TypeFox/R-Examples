##' @rdname show
##' @docType methods
##' @import methods
##' @aliases as.character
##' @param digits number of digits handed over to \code{format}
##' @param range should the values be indicated as range rather then first and
##'   last elements?
##' @param max.print maximum number of elements to be printed (of a variable)
##' @param shorten.to if a vector is longer than \code{max.print}, only the
##'   first \code{shorten.to[1]} and the last \code{shorten.to[2]} elements are
##'   printed
##' @return \code{as.character} returns a character vector fit to be printed by
##'   \code{cat} with \code{sep = "\n"}.
##' 
##' @seealso \code{\link[base]{as.character}}
##' @include paste.row.R
##' @export

setMethod ("as.character", signature = signature (x = "hyperSpec"),
           function (x, digits = getOption ("digits"), range = TRUE,
                     max.print = 5, shorten.to = c(2,1)){
  ## input checking
  validObject (x)

  if (is.null (max.print))
    max.print <- getOption ("max.print")

  if ((length (max.print) != 1) | ! is.numeric (max.print))
    stop ("max.print needs to be a number")
  if ((length (shorten.to) < 1) |(length (shorten.to) > 2) | ! is.numeric (shorten.to))
    stop ("shorten.to needs to be a numeric vector with length 1 or 2")
  if (sum (shorten.to) > max.print)
    stop ("sum (shorten.to) > max.print: this does not make sense.")

  ## printing information
  chr <- c("hyperSpec object",
           paste ("  ", nrow (x), "spectra"),
           paste ("  ", ncol (x), "data columns"),
           paste ("  ", nwl (x), "data points / spectrum")
           )

  chr <- c (chr, .paste.row (x@wavelength, x@label$.wavelength, "wavelength",
                             ins = 0, val = TRUE, range = FALSE,
                             shorten.to = shorten.to, max.print = max.print))
            
  n.cols <- ncol (x@data)

  chr <- c(chr, paste ("data: ", " (", nrow(x@data), " rows x ", n.cols,
                       " columns)", sep = ""))

  if (n.cols > 0)
    for (n in names (x@data))
      chr <- c(chr, .paste.row (x@data[[n]], x@label[[n]], n, ins = 3,
                                i = match (n, names (x@data)),
                                val = TRUE, range = range,
                                shorten.to = shorten.to, max.print = max.print))
  
  chr
})

##' Convert a hyerSpec object to character strings for Display
##' \code{print}, \code{show}, and \code{summary} show the result of
##' \code{as.character}.
##' 
##' \code{print}, \code{show}, and \code{summary} differ only in the defaults.
##' \code{show} displays the range of values instead,
##' @name show
##' @rdname show
##' @aliases show show,hyperSpec-method
##' @param object a \code{hyperSpec} object
##' @seealso \code{\link[methods]{show}}
##' @keywords methods print
##' @export
##' @examples
##' 
##' chondro
##' 
##' show (chondro)
##' 
##' summary (chondro)
##' 
##' print (chondro, range = TRUE)
##' 
setMethod ("show", signature = signature (object = "hyperSpec"), function (object){
  print (object, range = TRUE)
  invisible (NULL)
})

##'
##' \code{print} shows the overview giving the first and last values of each
##' data column (fastest).
##' @aliases print print,hyperSpec-method
##' @param x a \code{hyperSpec} object
##' @param ... \code{print} and \code{summary}  hand further arguments to \code{as.character}
##' @return \code{print} invisibly returns \code{x} after printing, \code{show} returns
##'   an invisible \code{NULL}.
##' @rdname show
##' @export 
##' @seealso \code{\link[base]{print}}
setMethod ("print", signature = signature (x = "hyperSpec"), function (x, range = FALSE, ...){
  validObject (x)
  cat (as.character (x, range = FALSE, ...), sep ="\n")
  invisible (x)
})


##' 
##' \code{summary} displays the logbook in addition.
##' 
##' @aliases summary summary,hyperSpec-method
##' @seealso \code{\link[base]{summary}}
##' @export
##' @rdname show
setMethod ("summary", signature = signature (object = "hyperSpec"),
           function (object, ...){
  print (object, ...)
})
