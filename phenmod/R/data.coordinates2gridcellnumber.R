data.coordinates2gridcellnumber <- function(grid, x,y){		
		valid.values <- which(is.na(grid$alt)==FALSE)
		valid.x <- grid$x[valid.values]
		valid.y <- grid$y[valid.values]

		test.x <- abs(valid.x - rep(x,length(valid.x)))
		lowest.value <- sort(test.x)[1]
		tocheck <- which(test.x == lowest.value)
		
		test.y <- abs(valid.y[tocheck] - rep(y, length(tocheck)))
		lowest.value <- sort(test.y)[1]
		found.y <- which(test.y == lowest.value)
	
		return(valid.values[tocheck[found.y]])
}