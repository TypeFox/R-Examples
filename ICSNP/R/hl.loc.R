hl.loc <- function(x, na.action=na.fail)
    {
    if (!is.vector(x)) stop("'x' must be a numeric vector")
    x<-na.action(x)
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    sums<-pair.sum(as.matrix(x))
    y.all<-c(x,sums/2)
    res <- median(y.all)
    return(res)
    }
