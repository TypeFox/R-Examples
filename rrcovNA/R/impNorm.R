impNorm <- function(x)
{
    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
                    dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dx <- dim(x)
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]

    ## impute the missing data using package norm
    s <- prelim.norm(x)                     # do preliminary manipulations
    thetahat <- em.norm(s, showits=FALSE)   # find the mle
    rngseed(1234567)                        # set random number generator seed
    ximp <- imp.norm(s, thetahat, x)        # impute missing data under the MLE
    xx<-imp.norm(s, thetahat, x)            # impute missing data under the MLE

    return(ximp)
}
