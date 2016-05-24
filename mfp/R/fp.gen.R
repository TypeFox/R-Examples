fp.gen <- function(x, pwrs, shift = NULL, scale = NULL)
{
    nobs <- length(x)
    pwr1 <- pwrs[1]
    pwr2 <- pwrs[2]
    X <- matrix(0, nrow = nobs, ncol = 2)   
    if(is.null(scale) | is.null(shift)) {
		x.transform <- fp.scale(x, scaling=TRUE)
		shift <- x.transform$shift
		scale <- x.transform$scale
		cat("[Pre-transformation : x->(x+", shift, ")/", scale, 
                "]\n", sep = "")
	}
	x <- x + shift
	x <- x/scale
#
# Deal with the first power
#
    if(!is.na(pwr1)) {
        x1 <- ifelse(pwr1 != rep(0, nobs), x^pwr1, log(x))
        X[, 1] <- x1
    }
#
# Other power
#
    if(!is.na(pwr2)) {
        if(pwr2 == pwr1)
            x2 <- log(x) * x1
        else x2 <- ifelse(pwr2 != rep(0, nobs), x^pwr2, log(x))
        X[, 2] <- x2
    }
    return(X)
}
