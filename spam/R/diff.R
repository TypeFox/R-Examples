# This is file ../spam/R/diff.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


########################################################################
diff.spam <- 
function (x, lag = 1, differences = 1, ...) 
{
    xlen <-   dim(x)[1L]
    if (length(lag) > 1L || length(differences) > 1L || lag < 
        1L || differences < 1L) 
        stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen) 
        return( numeric(0))

    for (i in 1L:differences){
      x <- x[(1L+lag):xlen,, drop = FALSE] - x[1L:(xlen-lag),, drop = FALSE] 
      xlen <- xlen - lag
    }
    return( x)
}

setMethod("diff","spam",diff.spam)
