dispdist <-
function(distout,ate,quantiles=c(0.025,0.975),display.plot=TRUE) {
	
	Pse <- function(x) mean((x - mean(x))^2)^.5
	
	greater.p.value <- mean(distout >= ate)
	lesser.p.value <- mean(distout <= ate)
	p.value <- 2*min(greater.p.value,lesser.p.value)
	p.value.alt <-  mean(abs(distout) >= abs(ate))
	
	quantileE <- quantile(distout,quantiles)
	sdE <- Pse(distout)
	exp.val <- mean(distout)
	
	if (display.plot == TRUE) {
		hist(distout,freq=TRUE,xlab="Estimated ATE",main=paste("Distribution of the Estimated ATE"),breaks=length(distout)^.5,lwd=1)
		lines(c(ate,ate),c(-1000,1000),col="red",lwd=2,lty=2)
		}
	return(list(two.tailed.p.value=p.value,two.tailed.p.value.abs=p.value.alt,greater.p.value=greater.p.value,lesser.p.value=lesser.p.value,quantile=quantileE,sd=sdE,exp.val=exp.val))
	}
