### bootstrap function for 2 by 2 contingency tables.
table.stats.boot <- function(CT, R = 100, alpha = 0.05, fudge = 0.01){
	OUT <- as.data.frame(matrix(NA, nrow = R, ncol = 4)  )
	names(OUT)<- c("pod", "far", "bias", "ets")

for(i in 1:R){
	N <- sum(CT) ## number of cases
CT.prob <- CT/N

L <- prod(dim(CT))  ### length of vector

X <- sample(1:L, size = N, replace = TRUE, prob = as.numeric(CT.prob))
   ### beware of zero entries.
CT.sample   <- matrix(tabulate(X, nbins = L),nrow = 2)
	temp    <- table.stats(CT.sample, silent = TRUE, fudge = fudge)
OUT$bias[i] <- temp$BIAS
OUT$pod[i]  <- temp$POD
OUT$far[i]  <- temp$FAR
OUT$ets[i]  <- temp$ETS	
	}
unOUT <- c(unlist(OUT))
if(any(is.nan(unOUT)) || any(is.na(unOUT)) || any(!is.finite(unOUT))) {
    wmsg <- paste("table.stats.boot: NaN, NA or non-finite numbers in one or more statistics.",
	           "Removing these values in calculating CIs.", sep="\n")
    warning(wmsg)
}
up <- apply(OUT,2,quantile, 1-alpha/2, na.rm=TRUE)
dw <- apply(OUT,2,quantile,  alpha/2, na.rm=TRUE)
return(rbind(up, dw))	
	}
	
