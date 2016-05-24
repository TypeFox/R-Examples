## -----------------------------------------------------------------------------
## Fonction Uniform
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

Uniform = function(x,y=NA,width=2) {
	if(is.na(y[1])){
		runif(length(x),min=-width/2,max=width/2)+x
	}
	else{
		1/width^length(x)*(max(abs(x-y))<=2)
	}
}
