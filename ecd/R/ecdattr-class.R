#' An S4 class to represent the \code{ecdattr} row in the Elliptic Database (ECDB)
#' 
#' The \code{ecdattr} class serves as an object-oriented interface 
#' between R and ECDB. This class is used extensively during
#' the \code{\link{bootstrap}} process. A list of light-weight \code{ecdattr} 
#' objects is created first by \code{\link{ecdattr.pairs}} function,
#' then the \code{\link{ecdattr.enrich}} function is invoked in parallel to calculate 
#' additional \code{ecd} attributes. 
#'
#' @slot call the match.call slot
#' @slot alpha numeric
#' @slot gamma numeric. When cusp is 1, \code{gamma} is derived.
#' @slot cusp  numeric, representing type of cusp. Only 0 (default) and 1 are allowed.
#' @slot use.mpfr logical, whether to use mpfr for ecd object.
#' @slot enriched logical. If \code{TRUE}, it indicates the object has been enriched with ecd attributes.
#' @slot alpha_m numeric, \code{alpha}*1000000.
#' @slot gamma_m numeric, \code{gamma}*1000000.
#' @slot ecd     an object of ecd class.
#' @slot attr    list of attributes. They are \code{NULL} protected for SQLite.
#'
#'
#' @include ecd-class.R
#' @keywords ecdattr
#'
#' @exportClass ecdattr
setClass("ecdattr",
         representation(call = "call",
                        alpha = "numeric",
                        gamma = "numeric",
                        cusp = "numeric",
                        use.mpfr = "logical",
                        enriched = "logical",
                        ecd = "ecd",
                        alpha_m = "numeric",
                        gamma_m = "numeric",
                        attr = "list"),
          prototype(call = call("ecdattr"),
                    alpha = NaN,
                    gamma = NaN,
                    cusp = 0,
                    use.mpfr = FALSE,
                    enriched = FALSE,
                    ecd = NULL,
                    alpha_m = NaN,
                    gamma_m = NaN,
                    attr = list())
)
### <---------------------------------------------------------------------->
