
# laser script to fitdAICrc.batch



fitdAICrc.batch <- function(x, observed=NULL, ...){
	rc1 <- apply(x, 1, pureBirth);
	rc2 <- apply(x, 1, bd);
	rv1 <- apply(x, 1, DDX);
	rv2 <- apply(x, 1, DDL);
	rv3 <- apply(x, 1, yule2rate, ...);	
	
	res <- matrix(0, nrow=nrow(x), ncol=6);
	for (i in 1:nrow(x)){
		res[i, 1] <- rc1[[i]]$aic;
		res[i, 2] <- rc2[[i]]$aic;
		res[i, 3] <- rv1[[i]]$aic;
		res[i, 4] <- rv2[[i]]$aic;
		res[i, 5] <- rv3[[i]]["AIC"];
		res[i, 6] <- min(rc1[[i]]$aic, rc2[[i]]$aic) - min(rv1[[i]]$aic, rv2[[i]]$aic, rv3[[i]]["AIC"]);
		
	}
	
	colnames(res) <- c('pb', 'bd', 'DDX', 'DDL', 'y2r', 'dAIC');
	res <- as.data.frame(res);
	if (!is.null(observed)){
		p.value <- length(res$dAIC[res$dAIC >= observed])/(nrow(res) + 1);
		cat('observed dAIC:', observed, '\n', 'p = ', p.value, '\n');
	}
	return(res);
}





















