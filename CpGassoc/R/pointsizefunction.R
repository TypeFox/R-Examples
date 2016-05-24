pointsizefunction <-
function(x) {
	y=1-x
	y[which(y<.2)]=.2
	y
	}
