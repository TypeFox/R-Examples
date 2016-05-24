# There are 25 numbered plot symbols; pch = plot character
xyplot( mcs ~ age, data = HELPrct, groups = sex, 
	    pch = c(1, 2), col = c('brown', 'darkgreen'), cex = .75 )  

