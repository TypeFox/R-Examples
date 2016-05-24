options(asciiType = "asciidoc")

##' Print ascii object
##'
##' Function displaying the asciidoc, txt2tags, reStructuredText, org or
##' textile code associated with the supplied object of class \code{ascii}.
##'
##' The package provides the new global option \code{asciiType}. Default value
##' is \code{"asciidoc"} (see examples).
##'
##' @param x An object of class \code{"asciiTable"}, \code{"asciiList"}, \code{"asciiMixed"}, \code{"asciiCbind"} or \code{"Report"}.
##' @param type Type of syntax produce.  Possible values for \code{type} are
##'   \code{"asciidoc"}, \code{"t2t"}, \code{"rest"}, \code{"org"},
##'   \code{"textile"} or \code{"pandoc"}.  Default value produce asciidoc syntax.
##' @param file A character string naming the file to print to. Default is
##'   \code{NULL} (print to the console).
##' @param append If \code{TRUE}, code will be appended to \code{file} instead
##'   of overwriting it. Default value is \code{FALSE}
##' @param escape If \code{TRUE}, characters in \code{list.escape} will be be
##'   printed with a \code{\\}. Default value is \code{FALSE}
##' @param list.escape Character vector. Default value is \code{c("\\\\_",
##'   "\\\\^")}
##' @param help logical print help? (objects of class \code{"Report"})
##' @param ... Additional arguments.  (Currently ignored.)
##' @author David Hajage \email{dhajage@@gmail.com}
##' @seealso \code{\link{ascii}}
##' @keywords print
##' @rdname print-ascii
##' @export
##' @examples
##' data(esoph)
##' ascii(esoph[1:10,])
##' print(ascii(esoph[1:10,]), type = "t2t")
##' print(ascii(esoph[1:10,]), type = "rest")
##' print(ascii(esoph[1:10,]), type = "org")
##' print(ascii(esoph[1:10,]), type = "textile")
##' print(ascii(esoph[1:10,]), type = "pandoc")
##' options(asciiType = "rest")
##' ascii(esoph[1:10,])
##' options(asciiType = "asciidoc")
setMethod("print","asciiTable",
           function(x, type = getOption("asciiType"), file = NULL, append = FALSE, escape = FALSE, list.escape = c("\\_", "\\^"), ...) {
             if (type == "asciidoc") res <- capture.output(x$show.asciidoc())
             if (type == "rest") res <- capture.output(x$show.rest())
             if (type == "org") res <- capture.output(x$show.org())
             if (type == "t2t") res <- capture.output(x$show.t2t())
             if (type == "textile") res <- capture.output(x$show.textile())
             if (type == "pandoc") res <- capture.output(x$show.pandoc())

             if (escape) {
               for (i in list.escape)
                 res <- gsub(i, paste("\\", i, sep = ""), res)
             }

             if (is.null(file)) {
               cat(res, sep = "\n")
             }
             else {
               if (append) op <- "a" else op <- "w"
               f <- file(file, op)
               writeLines(res, f)
               close(f)
             }
             invisible(x)
           }
           )

##' Show method for ascii objects
##'
##' @param object ascii or Report object
##' @rdname print-ascii
##' @export
setMethod(show, "asciiTable",
           function(object) {
             print(object)
           }
          )

##' @rdname print-ascii
##' @export
setMethod(print, "asciiList",
           function(x, type = getOption("asciiType"), file = NULL, append = FALSE, escape = FALSE, list.escape = c("\\_", "\\^"), ...) {
             if (type == "asciidoc") res <- capture.output(x$show.asciidoc())
             if (type == "rest") res <- capture.output(x$show.rest())
             if (type == "org") res <- capture.output(x$show.org())
             if (type == "t2t") res <- capture.output(x$show.t2t())
             if (type == "textile") res <- capture.output(x$show.textile())
             if (type == "pandoc") res <- capture.output(x$show.pandoc())

             if (escape) {
               for (i in list.escape)
                 res <- gsub(i, paste("\\", i, sep = ""), res)
             }

             if (is.null(file)) {
               cat(res, sep = "\n")
             }
             else {
               if (append) op <- "a" else op <- "w"
               f <- file(file, op)
               writeLines(res, f)
               close(f)
             }
             invisible(x)
           }
           )

##' @rdname print-ascii
##' @export
setMethod(show, "asciiList",
           function(object) {
             print(object)
           }
          )

##' @rdname print-ascii
##' @export
setMethod(print, "asciiMixed",
           function(x, type = getOption("asciiType"), file = NULL, append = FALSE, escape = FALSE, list.escape = c("\\_", "\\^"), ...) {
             if (type == "asciidoc") res <- capture.output(x$show.asciidoc())
             if (type == "rest") res <- capture.output(x$show.rest())
             if (type == "org") res <- capture.output(x$show.org())
             if (type == "t2t") res <- capture.output(x$show.t2t())
             if (type == "textile") res <- capture.output(x$show.textile())
             if (type == "pandoc") res <- capture.output(x$show.pandoc())

             if (escape) {
               for (i in list.escape)
                 res <- gsub(i, paste("\\", i, sep = ""), res)
             }

             if (is.null(file)) {
               cat(res, sep = "\n")
             }
             else {
               if (append) op <- "a" else op <- "w"
               f <- file(file, op)
               writeLines(res, f)
               close(f)
             }
             invisible(x)
           }
           )

##' @rdname print-ascii
##' @export
setMethod(show, "asciiMixed",
           function(object) {
             print(object)
           }
          )

##' @rdname print-ascii
##' @export
setMethod(print, "Report",
          function(x, help = FALSE, ...) {
            if (help)
              x$show.Report(help = TRUE)
            else
              x$show.Report(help = FALSE)
          })

##' @rdname print-ascii
##' @export
setMethod(show, "Report",
          function(object) {
            print(object)
          })
