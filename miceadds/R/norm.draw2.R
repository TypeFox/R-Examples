


#************************************************************
# use penalty parameter as the default
# in .norm.draw function
# This function is used in mice.impute.tricube.pmm2
.norm.draw2 <- function (y, ry, x, ridge = 1e-05, ...) 
{
    xobs <- x[ry, ]
    yobs <- y[ry]
    xtx <- t(xobs) %*% xobs
#    pen <- ridge * diag(xtx)
    pen <- ridge * ( 1+0* diag(xtx) )	
    v <- solve(xtx + diag(pen))
    coef <- t(yobs %*% xobs %*% v)
    residuals <- yobs - xobs %*% coef
	df_r <- sum(ry) - ncol(x)
	df_r <- max( df_r , 2 )
    sigma.star <- sqrt(sum((residuals)^2)/ stats::rchisq(1, df_r ))
    beta.star <- coef + (t(chol((v + t(v))/2)) %*% stats::rnorm(ncol(x))) * 
        sigma.star
    parm <- list(coef, beta.star, sigma.star)
    names(parm) <- c("coef", "beta", "sigma")
    return(parm)
}
#************************************************************
