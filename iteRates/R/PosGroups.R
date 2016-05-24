PosGroups <-
function(subTax,size){
	#wrapper - gives possible groupings of taxa of a given size
	temp<-combinations(length(subTax),size,subTax)
	return(temp)
	}

