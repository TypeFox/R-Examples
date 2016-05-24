summary.lad <-
function(object,...)
{
	lrtest <- function(object)
	{
		Tests <- Dfs <- Pvals <-cnames <- vector(length=object$numdir);

		for (w in 1:object$numdir)
		{
			Tests[w] <- 2*(object$loglik[object$numdir+1]-object$loglik[w]);
			Dfs[w] <- object$numpar[object$numdir+1]-object$numpar[w];
			Pvals[w] <- 1-pchisq(Tests[w], df=Dfs[w]);
			cnames[w] <- paste(paste(w-1, "D vs >= ", sep=""), paste(w, "D", sep=""), sep="");
		}
		ans <- data.frame(cbind(Tests, Dfs, Pvals)); rownames(ans) <- cnames;
		colnames(ans) <- c("Stat", "df", "p.value");

		return(ans)
	}

	ics <- function(object)
	{
		mat <- data.frame(rbind(object$aic, object$bic));
		v1 <- c("aic", "bic"); 
		v2 <-paste("d=", 0:(object$numdir), sep="");
		dimnames(mat) = list(v1, v2);
		return(mat)
	}
	ans <- object; class(ans) <-"summarylad";
	if (identical(object$numdir.test, TRUE))
	{	
		ans$LRT <- lrtest(object);
		ans$IC <- ics(object);
		return(ans);
	}
	return(ans)
}
