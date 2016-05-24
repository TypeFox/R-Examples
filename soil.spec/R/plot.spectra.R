#' Plot Spectra using both wavenumbers and wavelengths
#'
#' @author Andrew Sila \email{asila@cgiar.org} and Tomislav Hengl \email{tom.hengl@wur.nl}

plot.Spectra <- function(x,...){

	if(nrow(x@ab)>0){
	k <- ncol(x@ab)
	wavenumber <- as.numeric(sapply(names(x@ab[,-1]), function(obj){gsub("[^0-9.]", "", obj)}))  
	x <- as.matrix(x@ab[,-1])
	nano <- round(1e7/wavenumber, 2)
	wx <- summary(wavenumber)[-4]
	nx <- round(1e7/wx,2)
	xaxp <- c(min(wx),max(wx),length(wx)-1) ## Get x-axis tick marks
	par(mfrow=c(1,1))
	plot(wavenumber, x[1,], pch="", frame = FALSE, xlab="", ylab="", xlim=c(max(wavenumber), min(wavenumber)), ylim=c(min(na.omit(x)), max(na.omit(x))), xaxp=xaxp, ...)
	axis(3, wx, nx, lwd=0, tick=TRUE, lwd.ticks=0.1, col.ticks="black")
		for(i in 1: nrow(x)){
			lines(wavenumber, x[i,], col=i)
		  }
	box()
	mtext(expression("Wavenumbers cm"^-1), side=1, line=2, cex=0.8)
	mtext("Wavelength (nm)", side=3, line=2, cex=0.8)
	mtext("Absorbance", side=2, line=2, cex=0.8)
		} else {
			stop("Empty object")
		}
	}
setMethod("plot", signature("Spectra"), plot.Spectra)