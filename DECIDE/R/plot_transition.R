`plot_transition` <-
function(dataset) {
	relimp <- relative.importance(dataset)
	parameters <- relimp$parameters
	no_classes <- relimp$no_classes
	x1 <- seq(from=-4, to=4, length.out=1001)
	y <- list()
	for (i in 1:no_classes) {
		y[[i]] <- dnorm(x1, mean=parameters$mu[i], sd=parameters$sigma[i])
		}
	par(mar=c(5, 4, 4, 4) + 0.1)
    	plot(x1, y[[1]], xlim=c(-4,4), ylim=c(0,1), type="l", lwd=2, xlab="Academic performance", ylab="Proportion of cases")
	lines(x1, plogis(x1, location=-parameters$alpha[1]/parameters$beta[1], scale=1/parameters$beta[1]), xlim=c(-4,4), ylim=c(0,1), type="l", lty=2, lwd=2)
	for (i in 2:no_classes) {
		lines(x1, y[[i]], xlim=c(-4,4), ylim=c(0,1), type="l", lwd=2, col=i)
		lines(x1, plogis(x1, location=-parameters$alpha[i]/parameters$beta[i], scale=1/parameters$beta[i]), xlim=c(-4,4), 
				ylim=c(0,1), type="l", lty=2, lwd=2, col=i)
		}
	axis(4)
#    legend(x="topleft", legend=c("salariat", "intermediate", "working 
#	class", "salariat", "intermediate", "working class"), 
#	lty=c(1,1,1,2,2,2), col=c("black", "red", "green", "black", "red", "green"), cex=0.8, merge=T)

	mtext("Probability of transition", side=4, line=3)	
	}

