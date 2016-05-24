setMethod("extractAIC",
    signature(fit = "mle"),
    function (fit, scale, k = 2, ...) 
    {
 	   res <- logLik(fit)
 	   edf <- attr(res, "df")
	    c(edf, -2 * res + k * edf)
    }

)
