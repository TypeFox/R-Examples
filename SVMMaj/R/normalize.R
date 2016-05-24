transformdata <- function(x, standardize='interval', spline.knots=0, spline.degree=1){ 

		#nominal values --> binary columns
		if(is.factor(x))
      return(structure(model.matrix(~x-1),values = levels(x)))
		#numeric values --> splines or normalized values
		if(is.numeric(x))
      if(spline.knots!=0 || spline.degree>1 ) 
						return(isb(x,spline.knots=spline.knots,spline.degree=spline.degree))
				else	
						return(structure(normalize(x,standardize),
						dim= c(length(x),1)) )
		
		#logical values --> binary values
		if(is.logical(x))
        return(structure(as.numeric(x),
					 type='logical',
					 dim=c(length(x),1)))
					 
		#other values --> try to coerce into numeric values
		return(structure(as.numeric(x),
			dim = c(length(x),1)))
}
#=============================
#=============================
#=============================

predict.transDat <- function(x,attrib=NULL,
		values=NULL,standardization=NULL, splineInterval=NULL, splineDegree=NULL) {
	if(is.list(attrib)){
		values 			= attrib$values
		standardization	= attrib$standardization
		splineInterval	= attrib$splineInterval
		splineDegree	= attrib$splineDegree
	}
	#TRANSFORM AN ARRAY WITH PRESPECIFIED SETTINGS, IN PARTICULAR FROM TRANSFORMDATA
	if(!is.null(splineInterval))
		return(isb(x,knots=splineInterval,spline.degree=splineDegree))
	else if(!is.null(standardization)){
		I <- normalize(x,standardize=standardization)
		attr(I,'dim')<- c(length(x),1)
		return(I)
}	else if(!is.null(values)){
		x<- factor(x,levels=values)
		I<- model.matrix(~x-1)
		attr(I,'values')<-values
		return(I) 
}	else {
		I <- as.numeric(x)
		attr(I,'dim') <- c(length(x),1)
		return(I)
	}
}
#=============================
#=============================
#=============================

normalize <- function(x, standardize = 'zscore')
{
    #===================================================================
    #STANDARDIZE THE ATTRIBUTE MATRIX X
    #-------------------------------------------------------------------        
    if(is.list(standardize)) {
        stand.a <- standardize$a
        stand.b <- standardize$b
    } else if(standardize=='zscore') {
        stand.a <- mean(x)
        stand.b <- sd(x)
        stand.b[stand.b==0] <- 1
    } else if(standardize=='interval') {
        stand.a <- min(x)
        stand.b <- max(x) - stand.a
    } else {
        stand.a <- 0  
        stand.b <- 1
    }
    if(is.na(stand.b) || stand.b==0){
    	stand.a <- mean(x)
    	stand.b <-  1
    }

    X <- (x-stand.a)/stand.b
    attr(X,'standardization') <- list(a=stand.a,b=stand.b)
    return(X)
}

