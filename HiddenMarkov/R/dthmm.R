"dthmm" <-
function (x, Pi, delta, distn, pm, pn = NULL, discrete = NULL,
          nonstat = TRUE)
{
    if (is.null(discrete)){
        if (distn=="beta" | distn=="exp" | distn=="gamma" | distn=="lnorm" |
            distn=="logis" | distn=="norm") discrete <- FALSE
        else if (distn=="binom" | distn=="pois") discrete <- TRUE
        else stop("parameter discrete must be used when applying user distributions")
    }
    y <- c(list(x=x, Pi=Pi, delta=delta, distn=distn, pm=pm,
                pn=pn, discrete=discrete, nonstat=nonstat))
    class(y) <- "dthmm"
    return(y)
}

