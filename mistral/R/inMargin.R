## -----------------------------------------------------------------------------
## Fonction inMargin
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

inMargin = function(data,type,percent,alpha) {

	if(type=="SVM"){
		if(missing(alpha)) {alpha=1}
		isMargin = (abs(data$mean)<alpha)
	}
	if(type=="Kriging"){
		if(missing(alpha)) {
			if(missing(percent)) {alpha=1.96}
			else {alpha=qnorm((1+percent)/2)}
		}
		high = (data$mean+alpha*data$sd)
		low = (data$mean-alpha*data$sd)
		isMargin = (high*low)<0
	}

	return(isMargin)

} 
