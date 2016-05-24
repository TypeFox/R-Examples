plot.twostageTE <-
function(x, ...){	
	if(!inherits(x,"twostageTE")){
		stop("Error:  Object is not of class twostageTE")
	}
	plot_gpava <- function (x, main = "PAVA Plot", xlab = "Predictor", ylab = "Response", 
	col = "lightblue",...)
	{
	#x ... object of class pava

	  o <- order(x$z)
	  xval <- x$z[o]
	  yval <- x$x[o]
	  xcum <-  c(xval[1] - mean(diff(xval)), xval) 
	  jumps <- ((1:length(yval))[!duplicated(yval)]-1)[-1]   #jumps of fitted step function
	  jumps <- c(1, jumps, length(xval))
	  
	  lines(xval, yval, col = col, lwd = 1, type <- "S")
	  points(xval[jumps], yval[jumps], col = col, pch = 13)
	#  grid()
	}
	pava1 <- gpava(z=x$X1, y=x$Y1)
	if (!is.na(x$L2)) {
		pava2 <- gpava(z=x$X2, y=x$Y2)
	}

	if (!is.na(x$L2)) {
		plot(x=x$X1,y=x$Y1, pch="1", cex=1.5, xlab="", ylab="", ylim=range(c(x$Y1,x$Y2)), col="grey80")
		abline(h=x$threshold, lty=3, lwd=1, col=2)
		points(x=x$X2,y=x$Y2, pch="2", cex=1.5, col="grey65")
		plot_gpava(pava2, col="blue")
	}
	else {
		plot(x=x$X1,y=x$Y1, pch="1", cex=1.5, xlab="", ylab="", col="grey80")	
		abline(h=x$threshold, lty=3, lwd=1, col=2)
		plot_gpava(pava1, col=1)
	}
	abline(v=x$L1, lty=2, lwd=2)
	abline(v=x$U1, lty=2, lwd=2)
	if (!is.na(x$L2)) {
		abline(v=x$L2, col="blue", lwd=2)
		abline(v=x$U2, col="blue", lwd=2)
	}
	points(x=x$estimate, y=x$threshold, col="blue", pch=4, cex=1.5)
	if (!is.na(x$L2)) {
		segments(x$estimate,min(c(x$Y1,x$Y2))-1,x$estimate, x$threshold , lwd=2, col="blue")
	}
	else {
		segments(x$estimate,min(x$Y1)-1,x$estimate, x$threshold, lwd=2, col="blue")
	}	
	mtext("Explanatory", side=1, line=2.5, cex=1.65)
	mtext("Response", side=2 , line=2, cex=1.65)

	if (!is.na(x$L2)) {
		legend("topleft", c("Estimate", "1st Stage CI", "2nd Stage CI", "2nd Stage Iso-Regression"), pch=c(4, NA, NA, 13), col=c("blue",1, "blue","blue"), lty=c(NA,2,1,1), lwd=c(NA,2,2,1), bg="white")#, bty='n'
	}
	else {
		legend("topleft", c("Estimate", "1st Stage CI", "1st Stage Iso-Regression"), pch=c(4, NA, 13), col=c("blue",1,1), lty=c(NA,2,1), lwd=c(NA,2,1), bg="white")#, bty='n'
	}	
}
