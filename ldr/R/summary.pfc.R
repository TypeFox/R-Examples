summary.pfc <-
function(object,...)
{
  "%^%"<-function(M, pow) 
  { 
    if (prod(dim(M)==list(1,1))) return( as.matrix(M^pow) )
    eigenM = eigen(M) 
    return(eigenM$vectors%*%diag(c(eigenM$values)^pow)%*%t(eigenM$vectors))  
  }
  
	lrtest <- function(object)
	{
		Tests <- Dfs <- Pvals <-cnames <- vector(length=object$numdir)
		for (w in 1:object$numdir)
		{
			Tests[w] <- 2*(object$loglik[object$numdir+1]-object$loglik[w])
			Dfs[w] <- object$numpar[object$numdir+1]-object$numpar[w]
			Pvals[w] <- 1-pchisq(Tests[w], df=Dfs[w])
			cnames[w] <- paste(paste(w-1, "D vs >= ", sep=""), paste(w, "D", sep=""), sep="")
		}
		ans <- data.frame(cbind(Tests, Dfs, Pvals)) 
    rownames(ans) <- cnames
		colnames(ans) <- c("Stat", "df", "p.value")
		return(ans)
	}

	ics <- function(object)
	{
		mat <- data.frame(rbind(object$aic, object$bic))
		v1 <- c("aic", "bic") 
		v2 <-paste("d=", 0:(object$numdir), sep="")
		dimnames(mat) = list(v1, v2)
		return(mat)
	}

	r2 <- function(object)
	{
		R2 <- cnames <- vector(length=object$numdir)

		for (w in 1:object$numdir)
		{
			if (identical(object$numdir.test, FALSE)) tempd <-data.frame(cbind(object$y-mean(object$y), 
						object$Xc%*%((object$Deltahat)%^%(-1))%*%matrix(object$Gammahat[, 1:w], ncol=w))) 
			if (identical(object$numdir.test, TRUE)) tempd <-data.frame(cbind(object$y-mean(object$y), 
						object$Xc%*%((object$Deltahat[[w]])%^%(-1))%*%object$Gammahat[[w]])) 
			colnames(tempd)[1] <- "y"
			R2[w] <- summary(lm(y~.-1, data=tempd))$r.squared
		}
		Rsq <- data.frame(rbind(object$evalues, R2))
		v1 <- c("Eigenvalues", "R^2(OLS|pfc)") 
		v2 <-paste("Dir", 1:object$numdir, sep="")
		dimnames(Rsq) = list(v1, v2)
		return(Rsq)
	}
	ans <- object; class(ans) <-"summarypfc"

	if (identical(object$numdir.test, TRUE))
	{	
		ans$LRT <- lrtest(object)
		ans$IC <- ics(object)
		if (is.numeric(object$y)) ans$Rsq <- r2(object)
		return(ans)
	}
	return(ans)
}
