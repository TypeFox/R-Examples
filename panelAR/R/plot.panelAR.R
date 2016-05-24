plot.panelAR <- function(x,legend=TRUE,rot.axis=c(0,0),...){
	if(class(x)!="panelAR"){
		stop("x must be of class 'panelAR'.")
	}
	obs.mat <- x$panelStructure$obs.mat
	N.time <- ncol(obs.mat)
	N.panel <- nrow(obs.mat)
	times <- colnames(obs.mat)
	units <- rownames(obs.mat)
	image(z=t(obs.mat), ylab="", xlab="Time", axes=FALSE, col=c("darkred","wheat"),...)
    axis(1, at = seq(from = 0, to = 1, 
            length = N.time), tck = 0, labels=FALSE, lwd = 0, las = 2)
	text(seq(from = 0, to = 1, 
            length = N.time), par("usr")[3], labels = times, srt = rot.axis[1], pos = 1, 
            xpd =TRUE, cex=par("cex.axis"))
    axis(2, labels = FALSE, at = seq(from = 0, 
            to = 1, length = N.panel), tck = 0, lwd = 0, 
            las = 1)
    text(par("usr")[1],seq(from = 0, to = 1, 
            length = N.panel), labels =units, pos = 2, xpd = TRUE, cex=par("cex.axis"), 
            srt=rot.axis[2])
	if (legend) {
           par(xpd = TRUE)
            legend(x = 0.95, y = 1.1, col = c("darkred","wheat"), bty = "n", 
               xjust = 1, legend = c("Missing", "Observed"), 
               fill = c("darkred","wheat"), horiz = TRUE)
        }	
	}