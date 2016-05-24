omni.ate <-
function(Y,Z,perms,invert=FALSE,quantiles=c(0.025,0.975)) {
	
	Pse <- function(x) mean((x - mean(x))^2)^.5
		
	prob <- genprob(perms)
	ate <- estate(Y,Z,prob=prob)
	Ys0 <- genouts(Y,Z,0)
	suppressWarnings(dist0 <- gendist(Ys0,perms))
	greater.p.value <- mean(dist0 >= ate)
	lesser.p.value <- mean(dist0 <= ate)
	p.value <- 2*min(greater.p.value,lesser.p.value)
	p.value.alt <-  mean(abs(dist0) >= abs(ate))
	se.null <- Pse(dist0)
	
	hist(dist0,freq=TRUE,xlab="Estimated ATE",main=paste("Distribution of the Estimated ATE assuming tau = 0"),breaks=length(dist0)^.5,lwd=1)
	lines(c(ate,ate),c(-1000,1000),col="red",lwd=2,lty=2)
	
	YsA <- genouts(Y,Z,ate)
	suppressWarnings(distA <- gendist(YsA,perms))
	conf.int <- quantile(distA,quantiles)
	se <- Pse(distA)

	if (invert == FALSE) return(list(ate=ate,greater.p.value=greater.p.value,lesser.p.value=lesser.p.value,p.value=p.value,p.value.alt=p.value.alt,se.null=se.null,conf.int=conf.int,se=se)) else {
		conf.intInv <- c(invert.ci(Y,Z,prob,perms,quantiles[1]),invert.ci(Y,Z,prob,perms,quantiles[2]))
		return(list(ate=ate,greater.p.value=greater.p.value,lesser.p.value=lesser.p.value,p.value=p.value,p.value.alt=p.value.alt,se.null=se.null,conf.int=conf.int,se=se,conf.intInv=conf.intInv))
		
		
	}
}
