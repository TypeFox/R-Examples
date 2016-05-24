".as.area" <- function(x)
{
    ## Verifications
    if (!inherits(x, "data.frame"))
        stop("x should be of class \"data.frame\"")
    if (ncol(x) != 3)
        stop("x should have three columns")
    ## ID is again transormed into a factor
    if (!is.factor(x[,1]))
        x[,1] <-factor(x[,1])
    ## The class
    class(x)<-c("area", "data.frame")
    return(x)
  }

