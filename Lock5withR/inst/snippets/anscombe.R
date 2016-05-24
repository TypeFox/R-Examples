ddd <- with(anscombe, data.frame( 
	x = c(x1, x2, x3, x4), 
	y = c(y1, y2, y3, y4), 
	set = rep(1:4, each = nrow(anscombe))) 
	)
xyplot(y ~ x | factor(set), data = ddd, type = c('p', 'r'), as.table = TRUE, lty = 2)

