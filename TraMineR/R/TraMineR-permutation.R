TraMineR.permutation <- function(data, R, statistic, ...) {
	if (is.vector(data)) {
		n <- length(data)
	}
	else {
		n<- nrow(data)
	}
	perms <- list()
	perms$R <- sum(R)
	perms$t0 <- statistic(data,1:n,...)
	ntest <- length(perms$t0)
	if (perms$R>1) {
		perms$t <- matrix(NA, perms$R,length(perms$t0))
		perms$t[1,] <- perms$t0
		for(i in 2:(perms$R)){
			perms$t[i,] <- statistic(data,sample(n), ...)
		}
		perms$pval <- numeric(length(perms$t0))
	
		for(i in 1:length(perms$pval)){
			perms$pval[i] <- (sum(perms$t[, i]>=perms$t0[i]))/perms$R
		}
	}
	else {
		perms$t <- NA
		perms$pval <- as.numeric(rep(NA, ntest))
	}
	class(perms) <- "TraMineRPermut"
	return(perms)
}

TraMineR.permutationweight <- function(data, R, statistic, samplesize, sampleprob, t0, ...) {
	if (is.vector(data)) {
		n <- length(data)
	}
	else {
		n<- nrow(data)
	}
	perms <- list()
	perms$R <- sum(R)
	allind <- 1:n
	perms$t0 <- t0
	ntest <- length(perms$t0)
	# spop <- 1:samplesize
	#print(data.frame(allind,sampleprob))
	if (perms$R>1) {
		#perms$t0boot <- matrix(NA, perms$R,length(perms$t0))
		perms$t <- matrix(NA, perms$R,length(perms$t0))
		perms$t[1,] <- perms$t0
		#perms$t0boot[1,] <- perms$t0
		for(i in 2:perms$R){
			pop <- sort.int(sample(allind, size=samplesize, replace=TRUE, prob=sampleprob), method="quick")
			pop2 <- sample(allind, size=samplesize, replace=TRUE, prob=sampleprob)
			#bt <- statistic(pop,pop2, ...)
			perms$t[i,] <- statistic(pop,pop2, ...)
			#perms$t0boot[i,] <- bt[-(1:ntest)]
		}
		perms$pval <- numeric(length(perms$t0))
	
		for(i in 1:length(perms$pval)){
			perms$pval[i] <- sum(perms$t[, i]>=perms$t0[i])/perms$R
		}
	}
	else {
		perms$t <- NA
		perms$pval <- as.numeric(rep(NA, ntest))
	}
	class(perms) <- c("TraMineRPermut")
	return(perms)
}

summary.TraMineRPermut <- function(object, ...){
	perms <- object
	resume <- data.frame(t0 = perms$t0, "p-value"=perms$pval)
	rownames(resume) <- paste("t[,",1:length(perms$t0),"]")
	return(resume)
}
hist.TraMineRPermut <- function(x, index=1, breaks="FD",main=paste("Distribution of test statistic number ", index), xlab=index, pvalue.limit=NULL, freq=FALSE, type="density", txtcex=1, ...) {
	ntest <- length(x$t0)
	if(index >ntest){
		stop("index should be in range 1, ", ntest)
	}
	if (x$R==0) {
		stop("Cannot plot permutation test (none detected)")
	}
	if(type=="hist") {
		hist(x$t[, index], main=main, xlab=xlab, breaks=breaks, freq=freq, ...)
		if (!is.null(pvalue.limit)) {
			testbootorder <- order(x$t[, index], decreasing=TRUE)
			abline(v=x$t[testbootorder[round(pvalue.limit*x$R)], index], col="blue")
		}
		abline(v=x$t0[index], col="red")
	}
	else if(type=="density"){
		dens <- density(x$t[,index])
		##hist(x$t[,index], xlim = range(dens$x), xlab = "x", ylab = "density", freq = FALSE)
		plot(dens, main=main, xlab=xlab, ...)
		nullhypregion <- dens$x>=x$t0[index]
		polygon(x=c(min(dens$x[nullhypregion]),dens$x[nullhypregion]),y=c(0,dens$y[nullhypregion]),  col="red")
		if(txtcex>0) {
			text(dens$x[nullhypregion][1], dens$y[nullhypregion][1],labels=paste("  Statistic:",  format(x$t0[index], digits =3), "(P-value:",x$pval[index],")"), adj=c(0,0), cex= txtcex)
		}
	}

}

print.TraMineRPermut <- function(x,...){
	cat("\nPermutation test\n\n")
	print(summary(x),...)
	cat("P-value approximated using", x$R, "permutations\n")
}