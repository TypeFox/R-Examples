##' Sweave wrappers
##'
##' @param file Name of Sweave source file.
##' @param driver Sweave driver
##' @param syntax Sweave syntax
##' @param encoding Encoding
##' @param ... Further arguments passed to the driver's setup function.
##' @author David Hajage \email{dhajage@@gmail.com}
##' @seealso \code{\link{Sweave}}
##' @rdname sweave-wrapper
##' @export
##' @import utils
##' @keywords IO file
Asciidoc <- Sweave
formals(Asciidoc) <-  alist(file=, driver=RweaveAsciidoc, syntax=SweaveSyntaxNoweb, encoding="", ...=)

##' @rdname sweave-wrapper
##' @export
##' @import utils
##' @keywords IO file
T2t <- Sweave
formals(T2t) <-  alist(file=, driver=RweaveT2t, syntax=SweaveSyntaxNoweb, encoding="", ...=)

##' @rdname sweave-wrapper
##' @export
##' @import utils
##' @keywords IO file
ReST <- Sweave
formals(ReST) <-  alist(file=, driver=RweaveReST, syntax=SweaveSyntaxNoweb, encoding="", ...=)

##' @rdname sweave-wrapper
##' @export
##' @import utils
##' @keywords IO file
Org <- Sweave
formals(Org) <-  alist(file=, driver=RweaveOrg, syntax=SweaveSyntaxNoweb, encoding="", ...=)

##' @rdname sweave-wrapper
##' @export
##' @import utils
##' @keywords IO file
Textile <- Sweave
formals(Textile) <-  alist(file=, driver=RweaveTextile, syntax=SweaveSyntaxNoweb, encoding="", ...=)

##' @rdname sweave-wrapper
##' @export
##' @import utils
##' @keywords IO file
Pandoc <- Sweave
formals(Pandoc) <-  alist(file=, driver=RweavePandoc, syntax=SweaveSyntaxNoweb, encoding="", ...=)
