"beta.ratio" <-
function(a,b,x,y){
	### beta(a,b)/beta(x,y)
	exp( lgamma(a) + lgamma(b) - lgamma(a+b)
	-lgamma(x) - lgamma(y) + lgamma(x+y) )
}

