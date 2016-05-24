"perarea" <- function(x)
{
    ## Verifications
    if (!inherits(x, "area"))
        stop("x should be of class \"area\"")

    uu <- split(x[,2:3], x[,1])

    ## The function foo computes the perimeter of a polygon
    foo <- function(x) {
        if (!all(x[1,]==x[nrow(x),]))
            x <- rbind(x,x[nrow(x),])
        x1 <- x[-1,]
        x2 <- x[-nrow(x),]
        di <- sum(sqrt(((x2[,1]-x1[,1])^2)+((x2[,2]-x1[,2])^2)))
        return(di)
    }

    ## output
    res <- unlist(lapply(uu, foo))
    names(res) <- names(uu)
    return(res)
}

