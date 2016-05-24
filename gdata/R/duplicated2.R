duplicated2 <- function(x, bothWays=TRUE, ...)
    {
        if(!bothWays) {
            return(duplicated(x, ...))
        } else if(bothWays) {
            return((duplicated(x, ...) | duplicated(x, fromLast=TRUE, ...)))
        }
    }
