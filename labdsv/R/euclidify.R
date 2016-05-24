euclidify <- function (x,upper=FALSE,diag=FALSE)
{
    if (class(x) != 'dist') 
        stop('The first argument mmust be an object of class dist')
    tmp <- .Fortran("euclid",x=as.matrix(x),as.integer(attr(x,"Size")),
                     PACKAGE='labdsv')
    tmp2 <- as.dist(tmp$x)
    attr(tmp2, "Labels") <- dimnames(x)[[1]]
    attr(tmp2, "Diag") <- diag
    attr(tmp2, "Upper") <- upper
    attr(tmp2, "method") <- paste("euclidify", attr(x, "method"))
    attr(tmp2, "call") <- match.call()
    tmp2
}
