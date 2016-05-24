kappaTrap <-
function(x, l){
	if(l > 0) {
		ifelse(abs(x/l) > 2, 0, 
			ifelse(abs(x/l) < 1, 1, 2 - abs(x/l)))
	} else {
		ifelse(x == 0, 1, 0)
	}		
}

kappaRect <-
function(x, l){
	if(l > 0) {
		ifelse(abs(x/l) > 1, 0, 1)
	} else {
		ifelse(x == 0, 1, 0)
	}
}

kappaParzen <-
function(x, l) {
	if(l > 0) {
		ifelse(abs(x/l) < .5, 1 - 6*abs(x/l)^2 + 6*abs(x/l)^3, 
			ifelse(abs(x/l) < 1, 2*(1-abs(x/l))^3, 0))
	} else {
		ifelse(x == 0, 1, 0)
	}
}

kappaInDf <-
function(x, l) {
	if(l > 0) {
		ifelse(abs(x/l) < .6, 1, ifelse(abs(x/l) > 2, 0, 
			exp(-1*exp(-1/(abs(x/l)-.6)^2)/(abs(x/l)-2)^2)))
	} else {
		ifelse(x == 0, 1, 0)
	}
}