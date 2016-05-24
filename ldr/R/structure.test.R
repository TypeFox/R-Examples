structure.test <-
function(object1, object2)
{
	if (!(object1$structure %in% c("iso", "aniso", "unstr"))) 
	{
		struct <- object1$structure; temp <- paste("structure", struct, sep=" ")
		cat(paste(temp, " is not supported yet.", sep=" "))
		return()
	}

	if (!(object2$structure %in% c("iso", "aniso", "unstr"))) 
	{
		struct <- object2$structure; temp <- paste("structure", struct, sep=" ")
		cat(paste(temp, " is not supported yet.", sep=" "))
		return()
	}

	if (object1$numdir != object2$numdir)
	{
		cat("Models with different number of directions") 
		return()
	}

	if (!identical(object1$numdir.test, object2$numdir.test)) 
	{
		cat("Objects have different 'numdir.test' values") 
		return()
	}

	if (object1$numpar[1] > object2$numpar[1])
	{
		temp <- object2 
		object2<-object1 
		object1<-temp
	}

	if (!object1$numdir.test)
	{
		statistic <- 2*(object2$loglik - object1$loglik)
		thedf <- object2$numpar - object1$numpar
		pvalue <- 1 - pchisq(statistic, df=thedf)
		aics <- drop(c(object1$aic, object2$aic))
		bics <- drop(c(object1$bic, object2$bic))		
		LRT <- data.frame(cbind(statistic, thedf, round(pvalue, digits=4)))
		colnames(LRT) <- c("Stat", "df", "p.value"); rownames(LRT)=""
		AIC <- data.frame(rbind(aics))
		colnames(AIC) <- c(object1$structure, object2$structure); rownames(AIC)=""
		BIC <- data.frame(rbind(bics));
		colnames(BIC) <- c(object1$structure, object2$structure); rownames(BIC)=""
	}
	else 
	{
		d<- object1$numdir;
		if (max(object1$numdir) != max(object2$numdir)) d <- min(max(object1$numdir), max(object2$numdir))
		if (d==0){ message("There is no test with d=0"); return()}
		statistic <- thedf <- pvalue <- rnames <- vector(length=d)
		aics <- bics <- matrix(NA, ncol=2, nrow=d)
		for (i in 1:d)
		{
			statistic[i] <- 2*(object2$loglik[i+1] - object1$loglik[i+1])
			thedf[i] <- object2$numpar[i+1] - object1$numpar[i+1]
			pvalue[i] <- 1 - pchisq(statistic[i], df=thedf[i])
			rnames[i] <- paste("d = ", i, sep="")
			aics[i,] <- drop(c(object1$aic[i+1], object2$aic[i+1]))
			bics[i,] <- drop(c(object1$bic[i+1], object2$bic[i+1]))
		}
		LRT <- data.frame(cbind(statistic, thedf, round(pvalue, digits=4)))
		rownames(LRT) <-rnames;	colnames(LRT) <- c("Stat", "df", "p.value")
		AIC <- data.frame(rbind(aics))
		rownames(AIC) <- rnames; colnames(AIC) <- c(object1$structure, object2$structure)
		BIC <- data.frame(rbind(bics))
		rownames(BIC) <- rnames; colnames(BIC) <- c(object1$structure, object2$structure)
	}
	ans <- list(AIC=AIC, BIC=BIC, LRT=LRT)
	class(ans) <-"structure"
	return(ans) 
}
