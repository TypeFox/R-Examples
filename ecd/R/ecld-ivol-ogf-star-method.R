#' Calculate implied volatility using star OGF and small sigma formula
#'
#' Calculate implied volatility using star OGF and small sigma formula.
#' Only the postive ki and Lc is supported.
#' SGED is not supported yet.
#'
#' @param object an object of ecld class
#' @param ki a numeric vector of log-strike
#' @param epsilon numeric, small asymptotic premium added to local regime
#' @param otype option type
#' @param order.local numeric, order of the hypergeometric series to be computed
#'              for local regime. Default is \code{Inf}, use the incomplete gamma.
#'              When it is \code{NaN}, \code{L*} value is suppressed.
#' @param order.global numeric, order of the hypergeometric series to be computed
#'              for global regime. Default is \code{Inf}, use the incomplete gamma.
#' @param ignore.mu logical, ignore \code{exp(mu)} on both sides, default is \code{FALSE}.
#'
#' @return The state price of option in star OGF terms.
#'         For \code{ecld.ivol_ogf_star}, it is \eqn{\sigma_1}.
#'
#' @keywords ogf
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#' ld <- ecld(sigma=0.001)
#' ecld.ivol_ogf_star(ld, 0)
### <======================================================================>
"ecld.ivol_ogf_star" <- function(object, ki, epsilon=0, otype="c",
                        order.local=Inf, order.global=Inf, ignore.mu=FALSE)
{
    ecld.validate(object)

    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }

    if (is.na(object@mu_D)) object@mu_D <- ecld.mu_D(object)
    
    if (length(ki) > 1) {
        f <- function(ki) ecld.ivol_ogf_star(object, ki, epsilon, otype,
                    order.local=order.local,
                    order.global=order.global,
                    ignore.mu=ignore.mu)
        s <- ecld.mclapply(object, ki, f)
        return(s)
    }
    
    # ki should be length-one numeric
    stopifnot(length(ki)==1)
    
    one <- if(object@use.mpfr) ecd.mp1 else 1 # for gamma function
    
    lambda <- object@lambda * one
    sigma <- object@sigma * one
    
    k <- ki*object@sigma + object@mu
    R <- NaN 
    
    BS <- function(R) {
        sigma1 <- R*sigma
        mu_D1 <- -sigma1^2/4
        ki1 <- (k-mu_D1)/sigma1
        ld1 <- ecld(lambda=1, sigma=sigma1, mu=mu_D1)
        
        delta <- if (otype=="c" & ki<0) 1 else if (otype=="p" & ki>=0) -1 else 0
        dM1 <- delta * (exp(object@mu - object@mu_D)-1) * one + epsilon
        
        exp1 <- if (ignore.mu) 1 else exp(mu_D1)
        expL <- if (ignore.mu) 1 else exp(object@mu)
        
        # global regime
        Lc1 <- NULL
        if (is.na(order.global)) {
            stop("order.global can not be NA")
        } else if (order.global==Inf) {
            Lc1 <- sigma1 * exp1 * ecld.ogf_star(ld1, ki1)
        } else {
            Lc1 <- sigma1 * exp1 * ecld.ogf_star_exp(ld1, ki1, order.global)
        }
        
        # local regime
        Lc <- NULL
        if (is.na(order.local)) {
            Lc <- dM1
        } else if (order.local==Inf) {
            Lc  <- sigma  * expL * ecld.ogf_star(object, ki) + dM1
        } else {
            Lc  <- sigma  * expL * ecld.ogf_star_exp(object, ki, order.local) + dM1
        }
        return(Lc1 - Lc)
    }
    
    lower = 1
    upper = 100
    chk <- BS(lower)*BS(upper)
    if (!is.na(chk) & chk < 0) {
        rt <- ecd.uniroot(BS, lower=lower, upper=upper, use.mpfr=TRUE)
        R <- ecd.mp2f(rt$root)
    }
    return(R*object@sigma*one)
}
### <---------------------------------------------------------------------->
