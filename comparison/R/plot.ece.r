#############################################################################################
# functions to conduct emprical cross entropy calculations
# these functions are a reworking of Danial Ramos' ECE functions for Matlab
# (c) David Lucy 22 June 2010
# does the plotting
#############################################################################################
setMethod("plot",
signature(x = "ece"),
function(x, ...)
{


	# get rid of these eventually
	par(oma=c(0,0,0,0))
	par(mar=c(5,5,1,1))
	par(las=1)

	# eventually make these user selectable
	null.colour <- "black"
	ece.colour <- "red"
	calibrated.colour <- "blue"
	colours.vec <- c(null.colour, ece.colour, calibrated.colour)


	legend.text <- c("null", "observed", "calibrated")


	# dereference all the bits from the ece object
	prior 		<- slot(x, "prior")
	ece.null	<- slot(x, "ece.null")
	ece		<- slot(x, "ece")
	ece.cal		<- slot(x, "ece.cal")
	
	x.ordinates <- log10(prior / (1 - prior))
	legend.x.position <- x.ordinates[1]
	
	x.axis.text <- expression(paste(log[10], Odds(theta)))
	y.axis.text <- "empirical cross entropy"

	max.y <- max(c(ece.null, ece, ece.cal))
	
	plot(x.ordinates, ece.null, col=null.colour, ylim=c(0,max.y), type="l", xlab="", ylab="")

	points(x.ordinates, ece, col=ece.colour, type="l")
	points(x.ordinates, ece.cal, col=calibrated.colour, type="l")

	mtext(text=x.axis.text, side=1, line=2.5, cex=1.2)
	mtext(text=y.axis.text, side=2, line=2.8, cex=1.2, las=0)

	legend(legend.x.position, max.y, legend=legend.text, col=colours.vec, lty=1, bty="n")


	abline(h=0, lty=3)
	abline(v=0, lty=3)

}
)
