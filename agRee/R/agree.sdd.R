#---------------------------------------------------------------------------
agree.sdd <- function(ratings, conf.level=0.95, NAaction=c("fail", "omit")){

    if(!is.matrix(ratings) || ncol(ratings) < 2 || nrow(ratings) < 2)
      stop("'ratings' has to be a matrix of at least two columns and two rows.")

    na <- match.arg(NAaction)
    ratings <- switch(na,
                      fail = na.fail(ratings),
                      omit = na.omit(ratings))
    if(!is.matrix(ratings) || ncol(ratings) < 2|| nrow(ratings) < 2)
      stop("'ratings' has to be a matrix of at least two columns and two rows after removing missing values.")

    
    alpha <- 1 - conf.level
    k <- ncol(ratings)		
    n <- nrow(ratings)
    bar.x <- rowMeans(ratings)
    ssw <- sum((ratings - bar.x)^2)
    msw <- ssw / (n*(k-1))

    point <- sqrt(2)*qnorm(1-alpha/2)*sqrt(msw)
    lbound <- point*sqrt(n*(k-1)/qchisq(1-alpha/2, n*(k-1)))
    ubound <- point*sqrt(n*(k-1)/qchisq(alpha/2, n*(k-1)))
    list(value=point, lbound=lbound, ubound=ubound)
}
