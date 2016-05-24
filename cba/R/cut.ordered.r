
# cutting of ordinal variables
#
# ceeboo 2005

cut.ordered <- function(x, breaks, ...) {
    if (is.logical(breaks)) {
       if (length(breaks) != nlevels(x))
          stop("levels of",paste(sQuote("x"),"and",sQuote("breaks"),
               "do not conform"))
       breaks <- which(breaks)
    }
    else
       breaks <- sort(unique(breaks))
    if (is.character(breaks)) 
       breaks <- pmatch(breaks, levels(x))
    else
       breaks <- match(breaks, 1:nlevels(x))
    if (any(is.na(breaks)))
       stop(paste(sQuote("breaks"),"invalid"))
    #
    breaks <- unique(c(breaks, nlevels(x)))
    levels(x) <- rep(levels(x)[breaks], diff(c(0,breaks))) 
    x <- as.ordered(x)
    x
}

###
