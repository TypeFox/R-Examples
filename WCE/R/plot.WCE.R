plot.WCE <- function(x, allres = FALSE, ...){
	best <- which.min(x$info.criterion)
	sub1 = names(x$knotsmat)[best]
	if (x$a == FALSE) {sub2 = 'BIC'} else {sub2 = 'AIC'}
	n.knots <- length(x$info.criterion)
	if (allres == TRUE){
		matplot(t(x$WCEmat), lty=1, type = 'l', ylab = 'weights', xlab = 'Time elapsed') 
		title(paste('Estimated weight functions\n Best fit indicated with solid circles'))
		matplot(t(x$WCEmat), pch = 1, add = TRUE)
		points(x$WCEmat[best,], pch = 16, col = best)
		leg <- paste(names(x$knotsmat)[1], ',', sub2, '=', round(x$info.criterion[1],2))
			for (i in 2:n.knots){
				leg <- c(leg, paste(names(x$knotsmat)[i], ',', sub2, '=', round(x$info.criterion[i],2)))}
		if (x$constrained == 'Left') {
			legend('topleft', legend = leg, col = 1:n.knots, lty = 1)} 
		if (x$constrained == 'Right') {
			legend('topright', legend = leg, col = 1:n.knots, lty = 1)} 
		if (x$constrained == FALSE) {
			legend('topright', legend = leg, col = 1:n.knots, lty = 1)}
			} else {
		zetitle <- paste("Best-fitting estimated weight function\n", 
					paste(sub1, ' (', sub2, '=', round(x$info.criterion[best],2), ')', sep =''))
		plot(x$WCEmat[best,], ylab = 'weights', xlab = 'Time elapsed', pch=16)
		lines(x$WCEmat[best,], ylab = 'weights', xlab = 'Time elapsed')		
		title(zetitle)
		lines(x$WCEmat[best,])} 	
}