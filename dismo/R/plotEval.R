# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date :  December 2009
# Version 0.1
# Licence GPL v3


if (!isGeneric("plot")) {
	setGeneric("plot", function(x,y,...)
		standardGeneric("plot"))
}	


setMethod("plot", signature(x='ModelEvaluation', y='character'), 
	function(x, y='ROC', ext='max', main='', col='red', ...) {
		if (y == 'ROC') {
			txt = paste('AUC=', round(x@auc,3))
			plot(x@FPR, x@TPR, xlim=c(0,1), ylim=c(0,1), xlab='False postive rate', ylab='True positive rate', col=col, main=txt, ...)
			lines(x@FPR, x@TPR, col=col)
			lines(rbind(c(0,0), c(1,1)), lwd=2, col='grey')
		} else if (y %in% slotNames(x)) {
			dat <- slot(x, y)
			if (length(dat) != length(x@t)) {
				stop('not a valid slot for plotting')
			}
			if (ext=='max') {
				txt = paste('max at:', round(x@t[which.max(dat)], 2))
			} else if (ext=='min') {
				txt = paste('min at:', round(x@t[which.min(dat)], 2))
			} else { 
				txt = ''
			}
			if (main=='') { main <- y }
			if (txt != '') {
				main = paste(main, '-', txt)
			}
			plot(x@t, dat, xlab='threshold', ylab=y, col=col, main=main, ...)
			lines(x@t, dat, col=col)
		} else {
			stop('unknown value for "y"')
		}
	}
)


setMethod('density', signature(x='ModelEvaluation'), 
	function(x, ...) {
		pr <- density(x@presence)
		ab <- density(x@absence, bw=pr$bw )
		yl = c(min(ab$y, pr$y), max(ab$y, pr$y))
		xl = c(min(x@t), max(x@t))
		plot(ab, main='', ylab=paste('Density. Bandwidth=',round(pr$bw,5),paste=''), xlab='predicted value', xlim=xl, ylim=yl, lwd=2, lty=2, col='blue', ...)
		lines(pr, col='red', lwd=2)
#		x1 <- xl[1]+(xl[2]-xl[1]) / 3
#		x2 <- xl[1]+ 2 * (xl[2]-xl[1]) / 3
#		y = yl[1] + 0.5 * (yl[2] - yl[1])
#		text(x1,y,'absence',col='blue')
#		text(x2,y,'presence',col='red')
	} 
)

