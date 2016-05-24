dsvdis <- function(x, index, weight = rep(1,ncol(x)), step = 0., diag=FALSE, upper=FALSE)
{
    choices <- c("steinhaus", "sorensen", "ochiai", "ruzicka", "bray/curtis", "roberts", "chisq")
    i <- pmatch(index, choices)
    if(is.na(i))
        stop(paste(index, "is not a valid index:", paste(choices, 
            collapse = ", ")))
    if (!is.loaded("dsvdis")) {
        dyn.load("labdsv")
    }
    taxa <- deparse(substitute(x)) 
    x <- as.matrix(x)
    y <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    dis <- .Fortran("dsvdis",
        as.double(x),
        as.double(weight),
        as.integer(nrow(x)),
        as.integer(ncol(x)),
        as.integer(i),
        out = as.double(y),
        as.double(step),
        PACKAGE='labdsv')
    tmp <- matrix(dis$out, nrow = nrow(x))
    tmp2 <- tmp[row(tmp)>col(tmp)]
    class(tmp2) <- 'dist'
    attr(tmp2, "Labels") <- dimnames(x)[[1]]
    attr(tmp2, "Diag") <- diag
    attr(tmp2, "Upper") <- upper
    attr(tmp2, "method") <- choices[i]
    attr(tmp2, "call") <- match.call()
    attr(tmp2, "Size") <- nrow(x)
    attr(tmp2, "taxa") <- taxa
    return(tmp2)
}
