#' Constructor of ecld class
#' 
#' Construct an ecld class by providing the required parameters.
#' The default is the standard symmetric cusp distribution.
#' The default also doesn't calculate any ecd extension.
#' \code{ecld.from} allows you to pass the parameters from an existing ecd object.
#' \code{ecld.validate} checks if an object is ecld class.
#'
#' @param lambda numeric, the lambda parameter. Must be positive. Default: 3.
#' @param sigma numeric, the scale parameter. Must be positive. Default: 1.
#' @param beta  numeric, the skewness parameter. Default: 0.
#' @param mu    numeric, the location parameter. Default: 0.
#' @param with.ecd logical, also calculate the ecd object, default is \code{FALSE}.
#' @param with.mu_D logical, also calculate the risk-neutral drift, default is \code{FALSE}.
#'                  If \code{TRUE}, this flag supercedes \code{with.ecd}.
#'                  Also \code{mu} must set to zero.
#' @param with.RN logical, also calculate the risk-neutral ecd object, default is \code{FALSE}.
#'                If \code{TRUE}, this flag supercedes \code{with.mu_D}.
#' @param object an object of ecld class
#' @param is.sged logical, if \code{TRUE}, interpret parameters as SGED.
#' @param verbose logical, display timing information, for debugging purpose, default is \code{FALSE}.
#' @param sged.allowed logical, used in \code{ecld.validate} to indicate if the function allows SGED.
#' @param sged.only logical, used in \code{ecld.validate} to indicate if the function is only for SGED.
#'
#' @return an object of ecld class
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld
#' @export ecld.from
#' @export ecld.validate
#'
#' @examples
#' ld <- ecld()
#' ld <- ecld(2, 0.01)

### <======================================================================>
"ecld" <- function(lambda = 3, sigma = 1, beta = 0, mu = 0,
                  with.ecd = FALSE, with.mu_D = FALSE,
                  with.RN = FALSE, is.sged = FALSE, verbose=FALSE)
{
    call <- match.call()

    if(!is.numericMpfr(lambda)){
      stop("Parameter 'lambda' must be numericMpfr!\n")
    }
    if(!is.numericMpfr(sigma)){
        stop("Parameter 'sigma' must be numericMpfr!\n")
    }
    if(!is.numericMpfr(beta)){
        stop("Parameter 'beta' must be numericMpfr!\n")
    }
    if(!is.numericMpfr(mu)){
        stop("Parameter 'mu' must be numericMpfr!\n")
    }

    # -------------
    if(sigma <= 0){
        stop("Parameter 'sigma' must be positive!\n")
    }
    if(lambda <= 0){
        stop("Parameter 'lambda' must be positive!\n")
    }

    sum <- lambda+sigma+beta+mu
    if(!(abs(sum)>=0 & abs(sum) != Inf)) {
        stop(paste("Parameters must be finite, known real numbers!",
             "lambda=", ecd.mp2f(lambda), "sigma=", ecd.mp2f(sigma),
             "beta=", ecd.mp2f(beta), "mu=", ecd.mp2f(mu)
            ))
    }
    use.mpfr <- ifelse(class(sum)=="mpfr", TRUE, FALSE)

    if (use.mpfr) {
        sigma <- ecd.mp1 * sigma
    }

    ld <- new("ecld", call = call,
               lambda = unname(lambda),
               sigma = unname(sigma),
               beta  = unname(beta),
               mu    = unname(mu),
               use.mpfr = unname(use.mpfr),
               is.sged = is.sged,
               ecd = new("ecd"),
               ecd_RN = new("ecd"),
               status = 1L
          )
    
    # -------------
    # SGED
    if (is.sged) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: SGED"))
        if (with.ecd | with.mu_D | with.RN) {
            stop("SGED: with.ecd, with.mu_D, with.RN not supported")
        }
        return(invisible(ld))
    }
    
    # -------------
    # ecd
    if (with.ecd | with.mu_D | with.RN) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: calc ecd"))
        ld@ecd <- ecd(lambda=ecd.mp2f(lambda), beta=ecd.mp2f(beta),
                      sigma=sigma, mu=mu, verbose=verbose)
        ld@status <- bitwOr(ld@status, 2L)
    }
    if (with.mu_D | with.RN) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: calc mu_D"))
        ld@mu_D <- -log(ecd.imgf(ld@ecd, verbose=verbose))
        ld@status <- bitwOr(ld@status, 4L)
    }
    if (with.RN) {
        if (verbose) print(paste(Sys.time(), "ecld constructor: calc ecd_RN"))
        ld@ecd_RN <- ecd(lambda=ld@ecd@lambda, beta=ld@ecd@beta,
                         sigma=ld@ecd@sigma, mu=ld@mu_D, verbose=verbose)
        ld@status <- bitwOr(ld@status, 8L)
    }

    if (verbose) print(paste(Sys.time(), "ecld constructor: done"))
   
    invisible(ld)
}
### <---------------------------------------------------------------------->
#' @rdname ecld
"ecld.from" <- function(object,
                        with.ecd = FALSE, with.mu_D = FALSE,
                        with.RN = FALSE, verbose=FALSE)
{
    if (class(object) != "ecd") {
        stop("Must come from an ecd object")
    }
    if (object@alpha != 0 | object@gamma != 0) {
        stop("Must come from an ecd object with lambda-only parametrization")
    }
    
    ecld(lambda = object@lambda, sigma = object@sigma, beta = object@beta, mu = object@mu,
         with.ecd = with.ecd, with.mu_D = with.mu_D,
         with.RN = with.RN, verbose=verbose)
}
### <---------------------------------------------------------------------->
#' @rdname ecld
"ecld.validate" <- function(object, sged.allowed=FALSE, sged.only=FALSE)
{
    if (class(object) != "ecld") {
        stop(paste("Object must be an ecld object, found:", class(object)))
    }
    if (sged.only) {
        if (!object@is.sged) { # ECLD
            stop("Object must be an SGED object")
        } else {
            return(TRUE)
        }
    }
    if (object@is.sged) {
        if (sged.allowed) return(TRUE)
        stop("SGED object is not allowed here")
    }
    return(TRUE)
}
### <---------------------------------------------------------------------->

