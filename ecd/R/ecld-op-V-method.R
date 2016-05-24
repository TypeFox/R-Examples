#' The O, V, U operators in option pricing model
#'
#' The O operator takes a vector of implied volatility \eqn{\sigma_1(k)}
#' and transforms them to a vector of option state prices.
#' The V operator takes a vector of option state prices and transforms
#' them to a vector of implied volatility \eqn{\sigma_1(k)}.
#' The U operator calculates the log-slope of the option prices.
#'
#' @param sigma1 a vector of implied volatility (without T)
#' @param L a vector of option state prices
#' @param k a numeric vector of log-strike
#' @param sd numeric, the stdev of the distribution.
#'           Instead, if an ecld or ecd object is provided,
#'           the stdev will be calculated from it.
#' @param n numeric, number of lags in \code{ecld.op_U_lag}.
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#' @param stop.on.na logical, to stop if fails to find solution.
#'                   Default is to use NaN and not stop.
#' @param use.mc logical, to use mclapply (default), or else just use for loop.
#'               For loop option is typically for debugging.
#'
#' @return a numeric vector
#'
#' @keywords ogf
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.op_O
#' @export ecld.op_V
#' @export ecld.op_U_lag
#'
### <======================================================================>
"ecld.op_V" <- function(L, k, otype="c", stop.on.na=FALSE, use.mc=TRUE)
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    
    use.mpfr <- ifelse(class(L)=="mpfr" | class(k)=="mpfr", TRUE, FALSE)

    len <- length(L)
    if (len != length(k)) {
        stop("Length of L and k must match!")
    }
    if (len > 1 & use.mc==TRUE) {
        f <- function(i) ecld.op_V(L[i], k[i], otype=otype, stop.on.na=stop.on.na)
        s1 <- parallel::mclapply(1:len, f)
        s1 <- if (use.mpfr) ecd.mpfr(s1) else simplify2array(s1)
        return(s1)
    }
    # handle length-one numeric if use.mc
    
    s1 <- L*NaN
    for (i in 1:length(s1)) {
        df <- function(s) ecld.op_O(s, k[i], otype=otype) - L[i]
        
        lower = 0.01*ecd.mp1
        while (lower > 0.0000001*L[i] & df(lower) > 0) {
            lower <- lower/10
        }
        if (df(lower) > 0 & stop.on.na) {
            stop(paste("Failed to find starting lower for i=", i, "L=", L[i], "at k=", k[i]))
        }
        
        upper = L[i]*100**ecd.mp1
        if (upper < lower) upper <- lower*10
        while ( upper < 100 & df(upper) < 0) {
            upper <- upper*10
        }
        if (df(upper) < 0 & stop.on.na) {
            stop(paste("Failed to find starting upper for i=", i, "L=", L[i], "at k=", k[i]))
        }
        if (upper < lower & stop.on.na) {
            stop(paste("Failed to find starting lower/upper:", lower, upper, "for i=", i, "L=", L[i], "at k=", k[i]))
        }
        # print(c(i*ecd.mp1, L[i], k[i], lower, upper)) # debug

        # we use two iterations of uniroot to improve the precision to 10^-5 no matter how small sigma1 is
        if (df(lower) * df(upper) < 0 & upper > lower) {
            rt <- ecd.uniroot(df, lower=lower, upper=upper, use.mpfr=use.mpfr)
            s_appr <- rt$root
            df2 <- function(sn) ecld.op_O(sn*s_appr, k[i], otype=otype) - L[i]
            rt2 <- ecd.uniroot(df2, lower=0.9, upper=1.1, use.mpfr=use.mpfr)
            s1[i] <- s_appr * rt2$root
        } else {
            s1[i] <- NaN
        }
        if (is.na(s1[i]) & stop.on.na) {
            stop(paste("Failed to find root:", lower, upper, "for i=", i, "L=", L[i], "at k=", k[i]))
        }

    }
    return(s1)

}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_V
"ecld.op_O" <- function(sigma1, k, otype="c")
{
    if (!(otype %in% c("c","p"))) {
        stop(paste("Unknown option type:", otype))
    }
    use.mpfr <- ifelse(class(sigma1)=="mpfr" | class(k)=="mpfr", TRUE, FALSE)
    
    # use MPFR for erf
    sigma1 <- sigma1 * ecd.mp1
    M1k <- 1-exp(k)
    p <- -1/2 * ecd.erf(k/sigma1 - sigma1/4)
    q <- exp(k)/2 * ecd.erf(k/sigma1 + sigma1/4)

    sgn <- if (otype=="c") 1 else -1
    L <- p + q + sgn*M1k/2
    
    if (use.mpfr) return(L) else return(ecd.mp2f(L))
    
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_V
"ecld.op_U_lag" <- function(L, k, sd, n=2)
{
    if (class(sd)=="ecld") sd <- ecld.sd(sd)
    else if (class(sd)=="ecd") sd <- ecd.sd(sd)
    
    n2 <- floor(n/2)
    logL <- log(L)
    dlogL <- ecd.lag(logL, -n2) - ecd.lag(logL, n-n2)
    dk <- ecd.lag(k,-n2) - ecd.lag(k, n-n2)
    slope <- dlogL/dk*sd
    ecd.mp2f(slope)
}

