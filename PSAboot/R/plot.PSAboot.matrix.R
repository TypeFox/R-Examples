#' Matrix Plot of Boostrapped Propensity Score Analysis
#' 
#' @param bm result from \code{\link{PSAboot}}.
#' @export
matrixplot <- function(bm) {
	tmp <- reshape2::dcast(bm$pooled.summary[,c('iter','method','estimate')], 
				iter ~ method, value.var='estimate')
	panel.hist <- function(x, ...) {
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5) )
		h <- hist(x, plot = FALSE)
		breaks <- h$breaks; nB <- length(breaks)
		y <- h$counts; y <- y/max(y)
		rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
	}
	panel.density <- function(x, ...) {
		pu <- par("usr")
		d <- density(x,na.rm=TRUE,...)
		par("usr" = c(pu[1:2], 0, max(d$y)*1.5))
		lines(d)
		par("usr" = c(pu[1:2], 0, 1))
		boxplot(x, at=0.5, boxwex=0.3, horizontal=TRUE, add=TRUE)
	}
	panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		r <- abs(cor(x, y, use='pairwise.complete.obs'))
		txt <- format(c(r, 0.123456789), digits = digits)[1]
		txt <- paste0(prefix, txt)
		if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
		text(0.5, 0.5, txt, cex = cex.cor)
	}
	pairs(tmp[,2:ncol(tmp)], 
		  panel=panel.smooth, 
		  diag.panel=panel.density, 
		  upper.panel=panel.cor)
}
