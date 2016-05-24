pPval <- function(q, u, theta, sigma2){
    res <- rep(NA, length(q))
    for (j in 1:length(q)){
	   res[j] <- try(integrate(dPval, lower = 0, upper = q[j], u, theta, sigma2)$value, silent = TRUE)
	   if(class(res[j]) == "try-error"){res[j] <- 0.5}
	}
    return(res)
}

