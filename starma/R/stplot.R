stplot <- function(acf, ci, call, ggplot=T) {

	if (ggplot) {

		tlag <- value <- NULL	# avoid a note in "R CMD check" due to ggplot2 behaviour
	
		macf <- data.frame(slag=rep(colnames(acf), each=nrow(acf)), 
					 tlag=1:nrow(acf),
					 value=as.vector(acf))

		p <- ggplot(macf)
		p <- p + geom_linerange(ymin=0, aes(x=tlag, ymax=value))
		p <- p + geom_abline(intercept=0, slope=0)
		p <- p + geom_abline(intercept=ci, slope=0, color="blue", linetype="dashed")
		p <- p + geom_abline(intercept=-ci, slope=0, color="blue", linetype="dashed")
		p <- p + theme_bw()
		p <- p + scale_x_continuous(breaks=pretty_breaks())
		p <- p + theme(panel.grid=element_blank())
		p <- p + facet_wrap(~slag, ncol=1)
		p <- p + ggtitle(call)
	
		plot(p)

	} else {

		.pardefault <- par()

		par(mfrow=c(ncol(acf), 1), oma=c(.5,1.2,1,1), mar=c(3,2.5,0,0.8), 
			mgp=c(1.5,0.6,0))

		for (slag in 1:ncol(acf)) {
			plot(1:nrow(acf), acf[, slag], type="h", 
				xlab="T-lag", ylab=paste("S-lag =", slag - 1),
				cex.axis=.8,
				ylim=c(min(-1.1*ci, min(acf[, slag])), 
					max(1.1*ci, max(acf[, slag]))),
				xlim=c(0, nrow(acf) + 1), xaxs="i")
			lines(c(0, nrow(acf) + 1), c(0, 0))
			lines(c(0, nrow(acf) + 1), c(ci, ci), lty="dashed", col="blue")
			lines(c(0, nrow(acf) + 1), -c(ci, ci), lty="dashed", col="blue")
		}
	
		title(main=paste("Series"), outer=T)

		par(.pardefault)

	}

}