###################################################
#		Get transformed data
###################################################
Trans.Data <- function(n, div=10, pow.file="")
{
	res <- Get.Order.Div.Overdisp(n, div=div)
	ord <- res$ord
	x <- res$mea.log
	
	if (pow.file != "")
	{
		mat <- cbind(x, ord)
		colnames(mat) <- c('mean.log.mu', 'one.over.theta')
		write.table(mat, file=pow.file, row.names=F, col.names=T, quote=F)
	}
	
	### use linear regression
#	ord.pred <- predict(lm(ord ~ x, data=data.frame(x=x, ord=ord)), data.frame(x=log(rowMeans(n))))
	
	### use natrual cubic spline
	ord.pred <- predict(lm(ord ~ ns(x, df=3), data=data.frame(x=x, ord=ord)), data.frame(x=log(rowMeans(n))))
	
	n.trans <- n ^ (1 / ord.pred)
	
	ord <- Search.Ord(n.trans)
	cat("The power of the transformed data is:", ord, fill=T)
	
	return(n.trans)
}

###################################################
#		Get the overdispersion coefficient
###################################################
Get.Overdisp <- function(n, ord)
{
	SMALL.VAL <- 1e-8
	
	### take the transformation
	n1 <- n ^ ord
	
	### estimate the cmeans
	res <- Est.Depth(n1)
	cmeans <- res$cmeans
	keep <- res$keep
	
	### calculate the overdispersion
	n0 <- rowSums(n1) %*% t(cmeans)
	overdisp <- sum(rowSums((n1 - n0) ^ 2 / (n0 + SMALL.VAL))[keep]) - 
							(nrow(n) - 1) * (ncol(n) - 1) / 2
	
	return(overdisp)
}

###################################################
#		Search for the power
###################################################
Search.Ord <- function(n, st=0.5, ed=5.5, nint=5, sround=4)
{
	best.ord <- 1
	
	for (i in 1 : sround)
	{
		sseq <- seq(from=st, to=ed, by=(ed-st)/nint)
		overdisp <- rep(0, nint+1)
		for (j in 1 : (nint+1))
		{
			overdisp[j] <- Get.Overdisp(n, 1/sseq[j])
		}
		
		flag <- F
		for (j in 1 : nint)
		{
			if (overdisp[j] >= 0 & overdisp[j+1] <= 0)
			{
				st <- sseq[j]
				ed <- sseq[j + 1]
				
				best.ord <- switch((overdisp[j] + overdisp[j+1] >= 0) + 1, st, ed)
				flag <- T
				break
			}
		}
		
		if (!flag)
		{
			cat("ERROR in searching the power: the initial power range is too small!", fill=T)
		}
	}
	
	return(best.ord)
}

###################################################
#		Compare them: four divisions
###################################################
Get.Order.Div.Overdisp <- function(n, div)
{
	ord <- rep(0, div)
	mea.log <- rep(0, div)

	n <- n[order(rowSums(n), decreasing=T), ]
	
	nr <- floor(nrow(n) / div)
	for (i in 1 : div)
	{
		n.div <- n[((i - 1) * nr + 1) : (i * nr), ]
		
		mea.log[i] <- mean(log(rowMeans(n.div)))
		
		ord[i] <- Search.Ord(n.div)
	}	
	
	return(list(ord=ord, mea.log=mea.log))
}
