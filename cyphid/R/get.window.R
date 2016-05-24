get.window <-
function(dataset){
	xbar <- round(mean(dataset),1)
	xsd <- round(sd(dataset),1)
	upper <- xbar + (2*xsd)
	lower <- xbar - (2*xsd)
	win <- round(upper - lower,0)
	return(win)
	}

