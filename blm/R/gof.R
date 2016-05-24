gof <- function(object){
	
	if(class(object)[1]=="lexpit"|class(object)[1]=="blm"){
		p <- predict(object)
		y <- object@y
		w <- object@weights
	}
	else{
		p <- expit(predict(object)) # LOGISTIC
		y <- object$y
		w <- object$prior.weights
	}	
		brks <- quantile(p, seq(0, 1, length=11))
		g <- cut(p, br = brks, include.lowest = TRUE)
		
		O.cases <- tapply(w*y,g,sum)
		O.controls <- tapply(w*(1-y),g,sum)
		
		E.cases <- tapply(p*w,g,sum)
		E.controls <- tapply((1-p)*w,g,sum)

		result <- list(
							cases=data.frame(O=O.cases,E=E.cases), 
							controls=data.frame(O=O.controls, E=E.controls))
							
		chisq = sum((O.cases-E.cases)^2/E.cases+(O.controls-E.controls)^2/E.controls)

	    if(all(w==1)){
	        P = 1 - pchisq(chisq, df=8)
	        list(table = result, X2 = chisq, p.value = P)
	     }
	     else{
	     	# ADJUSTED F, ARCHER 2007
	     	df.num <- 9
	     	df.denom <- length(y)-8
	     	chisq <- chisq*((length(y)-8)/(10*length(y)))
	     	P <- 1 - pf(chisq,df.num,df.denom)
	     	list(table = result, X2 = chisq, p.value = P)
	     }

}

EO <- function(object, index=NULL,level=.95){
	if(class(object)[1]!="lexpit"&class(object)[1]!="blm")
		stop("Object must be an instance of a blm or lexpit model.")
	
	if(is.null(index)) index <- rep("Overall", length(object@y))
	
	p <- predict(object)
	
	O <- tapply(object@y, index, sum)
	E <- tapply(p*object@weights, index, sum)
	
	 z <- qnorm(1-(1-level)/2)
	 lower <- E/O*exp(-z*sqrt(1/O))
     upper <- E/O*exp(+z*sqrt(1/O))
      
  	results <- data.frame(
  					E = E, 
  					O = O, 
  					EtoO = E/O,
  				    lowerCI = lower, 
  				    upperCI = upper)
  				    
  	if(nrow(results)>1)
 		row.names(results) <- levels(index)

results
}