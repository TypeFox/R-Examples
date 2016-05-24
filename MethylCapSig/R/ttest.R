ttest <-
function(X, Y = NULL){
  ## ONE-SAMPLE t-TEST WHEN Y IS NULL
R <- rep(0,2);
if (is.null(Y)){
  	if (sum(X^2) != 0){
  	  	tmodel.rw <- t.test(x = colSums(X), y = NULL);
  	  	R <- c(tmodel.rw$statistic, tmodel.rw$p.value);
  	} else if (sum(X^2) == 0){
  		R[1] <- 0;
  		R[2] <- 1;
	}
} else if (!is.null(Y)){
	  ## REGIONWISE TEST - PERFORMED WHEN ALL THE OBSERVED SIGNALS ARE NON-ZERO
	  if (sum(X^2) + sum(Y^2) != 0){
		  tmodel.rw <- t.test(x = colSums(X), y = colSums(Y));
		  R <- c(tmodel.rw$statistic, tmodel.rw$p.value);
	  } else if (sum(X^2) + sum(Y^2) == 0){
	  	  R[1] <- 0;
	  	  R[2] <- 1;
	  }
}
return(R);
}
