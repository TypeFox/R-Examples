"raldous" <-
function(tip.number){

	if (tip.number < 1 | tip.number!=floor(tip.number)){
		stop("tip.number must be an integer greater than 1")
	}

	res=rtreeshape(n=1, tip.number=tip.number, model="aldous")
	res

}

