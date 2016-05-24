#' Compute degrees of freedom for a 2-sample t-test
#' 
#' This function computes degrees of freedom for a 2-sample t-test from the standard deviations 
#' and sample sizes of the two samples.
#' 
#' @param sd1 standard deviation of the sample 1
#' @param sd2 standard deviation of the sample 2
#' @param n1 size of sample 1
#' @param n2 size of sample 2
#' @return estimated degrees of freedom for 2-sample t-test
#' @export
#' @examples
#' data(KidsFeet, package="mosaicData")
#' fs <- favstats( length ~ sex, data=KidsFeet ); fs
#' t.test( length ~ sex, data=KidsFeet )
#' tdf( fs[1,'sd'], fs[2,'sd'], fs[1,'n'], fs[2,'n'])
#' 
tdf <-
function (sd1, sd2, n1, n2) 
{
    v1 <- sd1^2/n1
    v2 <- sd2^2/n2
    df1 <- n1 - 1
    df2 <- n2 - 1
    num <- (v1 + v2)^2
    denom <- (v1^2/df1) + (v2^2/df2)
    return(num/denom)
}
