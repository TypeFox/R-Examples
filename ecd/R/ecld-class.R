#' An S4 class to represent the lambda distribution
#' 
#' The \code{ecld} class serves as an object-oriented interface for the lambda distribution. 
#' The \code{ecld} prefix will also be used as the namespace for many analytic formulai
#' derived in lambda distribution, especially when lambda = 1,2,3.
#' Because of the extensive use of analytic formulai and enhanced precision through
#' the unit distribution, MPFR is not needed in most cases. This makes option pricing
#' calculation in \code{ecld} much faster than its counterpart built on the more
#' general-purpose \code{ecd} library.
#'
#' @slot call the match.call slot
#' @slot lambda numeric
#' @slot sigma numeric
#' @slot beta  numeric
#' @slot mu  numeric
#' @slot use.mpfr logical, whether to use mpfr for ecld object. If any of the above parameters
#'                is mpfr, then this flag is set to \code{TRUE}.
#' @slot is.sged logical, if \code{TRUE}, interpret parameters as SGED.
#' @slot ecd     the companion object of ecd class (optional)
#' @slot mu_D    the risk-neutral drift (optional)
#' @slot ecd_RN  the risk-neutral companion object of ecd class (optional)
#' @slot status  numeric, bitmap recording the state of the calculation layers.
#'               1: bare bone; 2: ecd; 4: mu_D; 8: ecd_RN
#'
#' @include ecd-class.R
#' @keywords ecld
#' 
#' @author Stephen H. Lihn
#'
#' @exportClass ecld
setClass("ecld",
         representation(call = "call",
                        lambda = "numericMpfr",
                        sigma = "numericMpfr",
                        beta = "numericMpfr",
                        mu = "numericMpfr",
                        use.mpfr = "logical",
                        is.sged = "logical",
                        ecd = "ecd",
                        mu_D = "numericMpfr",
                        ecd_RN = "ecd",
                        status = "numeric"),
          prototype(call = call("ecld"),
                    lambda = NaN,
                    sigma = NaN,
                    beta = NaN,
                    mu = NaN,
                    use.mpfr = FALSE,
                    is.sged = FALSE,
                    ecd = NULL,
                    mu_D = NaN,
                    ecd_RN = NULL,
                    status = 0)
)
### <---------------------------------------------------------------------->
