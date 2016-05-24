qq4res <- function(res, main = '', xlabel = 'Theoretical Quantiles', ylabel = 'Residuals', lgd = NULL){
	qqnorm(res, main = main, ylab = ylabel, cex = 2.5, pch = 16, cex.lab = 2, cex.axis = 2) 
	qqline(res, cex = 3)
	
	if (is.null(lgd) == FALSE) {
		legend('topleft', lgd, cex = 2)
	}
}
