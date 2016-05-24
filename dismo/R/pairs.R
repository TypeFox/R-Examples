 
 if (!isGeneric("pairs")) {
	setGeneric("pairs", function(x, ...)
		standardGeneric("pairs"))
}
 

setMethod('pairs', signature(x='DistModel'), 
	function(x, v=NULL, pa='pa', hist=TRUE, cor=TRUE) {
	
		panelhist <- function(x,...)	{
			usr <- par("usr")
			on.exit(graphics::par(usr))
			graphics::par(usr = c(usr[1:2], 0, 1.5) )
			h <- hist(x, plot = FALSE)
			breaks <- h$breaks
			nB <- length(breaks)
			y <- h$counts
			y <- y/max(y)
			rect(breaks[-nB], 0, breaks[-1], y, col="green")
		}
		
		panelcor <- function(x, y,...) {
			usr <- par("usr")
			on.exit(graphics::par(usr))
			graphics::par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			txt <- format(c(r, 0.123456789), digits=2)[1]
			text(0.5, 0.5, txt, cex = max(0.5, r * 2))
		}
	
		if (hist) {dp <- panelhist} else {dp <- NULL}
		if (cor) {up <- panelcor} else {up <- NULL}
	
	
		if (is.null(v)) {
			v <- 1:length(colnames(x@presence))
		} 
		if (length(v) < 2) {
			stop('pairs needs at least 2 variables')
		}

		if (! x@hasabsence) { pa <- 'p' }
		
		padj <- 4/length(v)
		if (pa=='pa') {
			d <- rbind(x@absence, x@presence)[,v]
			i <- c(rep(1, nrow(x@absence)), rep(2, nrow(x@presence)))
			pairs(d, main='', pch = c(4,21)[i], cex=c(0.01, padj)[i], bg=c('black','red')[i], upper.panel=up, diag.panel=dp)
		} else {
			if (pa=='p') { 
				pairs(x@presence[,v], main='presence', pch = 21, cex=padj, bg='red', upper.panel=up, diag.panel=dp)
			} else {
				pairs(x@absence[,v], main='absence', pch = 21, cex=padj, bg='red', upper.panel=up, diag.panel=dp)
			}
		}
	}
)

