##' @rdname read.jdx
##' @export
##' @note \code{read.jdx.Shimadzu}  is now defunct. Please use \code{read.jdx} instead.
read.jdx.Shimadzu <- function (...){
  .Defunct ("read.jdx",
            package = "hyperSpec",
            msg = "read.jdx.Shimadzu is now defunct.\nPlease use\nread.jdx (filename, encoding, header = list (xunits = expression (m/z), yunits = 'I / a.u.'))\ninstead.")
}

