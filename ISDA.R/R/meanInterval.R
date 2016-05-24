meanInterval <-
function(interval){
	xmin = interval$minValue
	xmax = interval$maxValue
	m =length(xmin)
	mean = 0
	if(m != length(xmax)){
		print("the lenght of xmin and xmax must be the same")
	}else{
		mean = (sum(xmin)+ sum(xmax))/(2*m)
	}
	mean
}

