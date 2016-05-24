if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')

int2alpha <- function(i, alpha)
{
    ## {0,1,..} |--> "alphabet" representation of discrete t.s.
    ss(alpha)[1L + i]
}
alpha2int <- function(x, alpha)
{
    ## (single) character |--> {0,1,..}  representation of discrete t.s.
    match(x, ss(alpha)) - 1L
}

int2char <- function(i, alpha)
{
    ## {0,1,..} |--> "alphabet" representation -- as 1 string --
    paste(int2alpha(i,alpha), collapse="")
}
char2int <- function(x, alpha)
{
    ## 1-string |--> {0,1,..}  representation of discrete t.s.
    match(ss(x), ss(alpha)) - 1L
}

id2ctxt <- function(id, m = nchar(alpha), alpha = NULL) {
    ## Purpose: Compute context from "ID" as returned by predict.vlmc

    if((m <- as.integer(m)) < 2)
        stop("alphabet length 'm' must be >= 2")
    ## Improve (but then, use C anyway!):
    r <- vector("list", n <- length(id <- as.integer(id)))
    i.ok <- !is.na(id)
    r[!i.ok] <- NA
    lev <- floor(1e-7 + log(id, m))

    for(i in 1:n) if(i.ok[i]) {
        ii <- id[i]
        ## convert ID 'ii' to {0:(m-1)} coded vector 'rr':
        rr <- integer(lev[i])
        for(ll in lev[i]:1) {
            rr[ll] <- ii %% m
            ii <- ii %/% m
        }
        r[[i]] <- rr
    }
    if(is.null(alpha) || (is.logical(alpha) && !alpha))
        r # list of integer vectors
    else if(is.logical(alpha) && alpha)
        ## return string, using "01.." alphabet
        sapply(r, function(i)paste(i, collapse=""))
    else if(is.character(alpha)) { ## using  'alpha' alphabet
        if(length(alpha) > 1) ## return vector of characters
            sapply(r, function(i) alpha[1L + i])
        else ## return string
            sapply(r, function(i) int2char(i,alpha))
    }
    else {
        warning("invalid 'alpha'; using alpha = NULL")
        r
    }
}
## })# local
