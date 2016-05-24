#########################################################################
# PLOT PROJECTED POPULATION
#########################################################################

plot_projection <- function(projection, kind) {

	if (missing(kind)) {

		mydata <- projection$runs
		dyears <- 1:dim(mydata)[1]
		dval <- apply(mydata, 3, rowSums)
		dcl <- colorRampPalette(c("blue", "yellow"))(ncol(dval))
		plot(0, 0, type="l", xaxt="n", xlab="Time steps", ylab="Individuals", xlim=c(1,nrow(mydata)), ylim=c(0, max(dval)), main="Total population size")
		for (i in 1:ncol(dval)) {
			lines(dyears, dval[,i], col=dcl[i])
		}
		axis(1, at=dyears, labels=dyears, las=1)


	} else {

		mydata <- colSums(aperm(projection$runs, c(2,1,3)))
		dyears <- 1:nrow(mydata)

		if (kind == "mean") dmed <- apply(mydata, 1, mean)
		if (kind == "median") dmed <- apply(mydata, 1, median)

		dlci <- apply(mydata, 1, quantile, 0.025)
		duci <- apply(mydata, 1, quantile, 0.975)

		plot(0,0, type="l", xaxt="n", xlab="Time steps", ylab="Individuals", xlim=c(1,nrow(mydata)), ylim=c(0, max(duci)), main="Total population size")
		polygon(x=c(dyears, rev(dyears)), y=c(dlci, rev(duci)), col="gray95", border=NA)
		lines(dyears, dlci, lty="dashed")
		lines(dyears, duci, lty="dashed")
		lines(dyears, dmed, lty="solid")
		points(dyears, dmed[dyears], pch=21, bg="white", cex=1.5)
		axis(1, at=dyears, labels=dyears, las=1)

	}

}


