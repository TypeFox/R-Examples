formatPval <- function (pv,
                        digits = max(1, getOption("digits") - 2),
                        eps = 0.0001,
                        na.form = "NA",                                           
                        scientific=FALSE,
                        includeEquality=FALSE) 
{
    ## first discard NA values, which will be included as the string in "na.form"
    ## at the end of the function 
    if ((has.na <- any(ina <- is.na(pv))))
    {
        pv <- pv[!ina]
    }

    r <- character(length(is0 <- pv < eps))

    ## process the large p values
    if (any(! is0))
    {
        rr <- pv <- pv[! is0]

        expo <- floor(log10(ifelse(pv > 0, pv, 1e-50)))
        fixp <- expo >= -3 | (expo == -4 & digits > 1)

        if (any(fixp))
        {
            ## DSB's initial version:
            rr[fixp] <- format(pv[fixp], digits=digits, scientific=scientific)
            
            ## my version:
            rr[fixp] <- disp(pv[fixp], 2, 2)
        }

        if (any(!fixp))
        {
            ## DSB's initial version:
            rr[! fixp] <- format(pv[! fixp], digits=digits, scientific=scientific)

            ## my version:
            rr[! fixp] <- disp(pv[! fixp], 2, 2)
        }

        r[! is0] <- rr
    }

    ## process the small p values
    if (any(is0))
    {
        digits <- max(1, digits - 2)
        
        if (any(!is0))
        {
            nc <- max(nchar(rr, type = "w"))

            if (digits > 1 && digits + 6 > nc)
            {
                digits <- max(1, nc - 7)
            }
        }

        r[is0] <- format(eps, digits = digits, scientific=scientific)
    }

    ## add (in)equality signs
    frontEqual <- 
        if(includeEquality)
            "= "
        else
            ""
    r <- paste(ifelse(is0, "< ", frontEqual),
               r,
               sep="")
    
    ## finally add back the NAs
    if (has.na)
    {
        rok <- r
        r <- character(length(ina))
        r[! ina] <- rok
        r[ina] <- na.form
    }
    
    return(r)
}
