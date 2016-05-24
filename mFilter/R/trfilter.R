###
### Trigonometric regression filter
trfilter <- function(x,pl=NULL,pu=NULL,drift=FALSE)
{
    if(is.null(drift)) drift <- FALSE
    call <- as.call(match.call())
    xname <- deparse(substitute(x))
    res <- cffilter(x,pl=pl,pu=pu,drift=drift,root=FALSE,
                    type="trigonometric",nfix=0,theta=1)
    res$method <- "trfilter"
    res$call <- call
    res$xname <- xname
    return(res)
}

