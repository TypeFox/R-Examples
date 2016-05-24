"mmglmlong1" <-
function (y, Pi, delta, glmfamily, beta, Xdesign, longitude,
          sigma=NA, nonstat=TRUE, size=NA, msg=TRUE){
    if (msg)
       message("NOTE: 'mmglmlong1' and its methods are still under development and may change.")
    if (glmfamily$family=="binomial" & all(is.na(size)))
        stop("Argument size must be specified when fitting a binomial model.")
    x <- c(list(y=y, Pi=Pi, delta=delta, glmfamily=glmfamily,
                beta=beta, Xdesign=Xdesign, sigma=sigma,
                nonstat=nonstat, longitude=longitude, size=size))
    class(x) <- c("mmglmlong1")
    return(x)
}

