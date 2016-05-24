reduceVar <-
function(x, breaks=NULL){
    if(is.numeric(x) | is.integer(x)){
        if(is.null(breaks)){x
            breaks <- "sturges"
        }
        if(is.character(breaks)){
            breaks <- match.arg(tolower(breaks), c("sturges",
                                                   "fd", "scott", "ss"))
            breaks <- switch(breaks, sturges = nclass.Sturges(x),
                             fd = nclass.FD(x),
                             scott = nclass.scott(x),
                             ss = nclass.ss(x),
                             stop("unknown 'breaks' algorithm"))
        }
        if(length(breaks) > 0){
            if(length(breaks)==1){
                rg <- range(x, na.rm=TRUE)
                breaks <- seq(rg[1],rg[2], length = breaks)
            }
            breaks <- unique(breaks)
            if(length(breaks)>1)
                x <- cut(x, breaks=breaks, include.lowest = TRUE, labels = FALSE)
            else
                x <- as.numeric(x)
        }
    } else {
        x <- as.numeric(x)
    }
    return(list(x=x, breaks=breaks))
}
