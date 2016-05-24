`rangeCuts` <-
function(x,n = 5,params = NA) {
	labs = 1:(n-1)
	frac = (1:(n-1))/n
	res = min(x)+frac*(max(x)-min(x))
	names(res) = paste("#",labs,sep="")
	res}

