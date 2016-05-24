`covStruct.create` <- 
function(covtype, d, known.covparam, var.names, coef.cov=NULL, coef.var=NULL, nugget=NULL, nugget.estim=FALSE, nugget.flag=FALSE, iso=FALSE, scaling=FALSE, knots=NULL, kernel=NULL) {
	
  if (covtype=="matern5_2add0") {
    weight <- coef.cov[(d+1):(2*d)]
    covStruct <- new("covAdditive0", 
                     d = as.integer(d),  
                     name = "matern5_2add0",
                     var.names = as.character(var.names),
                     sd2 = as.numeric(sum(weight)),
                     known.covparam = as.character(known.covparam),
                     range.val = as.numeric(coef.cov[1:d]),
                     range.names = paste("range", var.names, sep="."),
                     weight = as.numeric(weight),
                     weight.names = paste("weight", var.names, sep="."),
                     nugget = as.numeric(nugget),
                     nugget.flag = TRUE,
                     nugget.estim = nugget.estim,
                     param.n = as.integer(2*d+1)
    )
    return(covStruct)
  }
  
  if( covtype=="covUser" ){
    covStruct <- new("covUser", kernel=kernel, nugget.flag=length(nugget)>0, nugget=as.double(nugget))
	  return(covStruct)
  }
  
  if (scaling & iso) {
    iso <- FALSE
    warning("At this stage no isotropic version is available, regular scaling is applied.")
  }
  
	covsetI <- c("gauss", "exp", "matern3_2", "matern5_2")
	covsetII <- c("powexp")

	classType <- "covTensorProduct"
	if (iso) classType <- "covIso"
	if (scaling) {
    if (is.null(knots)) {
      classType <- "covAffineScaling"
    } else classType <- "covScaling"
	}
  
	covStruct <- new(classType, d=as.integer(d), name=as.character(covtype), 
		sd2 = as.numeric(coef.var), var.names=as.character(var.names), 
		nugget = as.double(nugget), nugget.flag=nugget.flag, nugget.estim=nugget.estim, known.covparam=known.covparam) 
									                
	if (!scaling) {					
		
		covStruct@range.names  = "theta"
		
		if (is.element(covtype, covsetI)) {
			covStruct@paramset.n <- as.integer(1)
			if (iso) {
				covStruct@param.n <- as.integer(1)
			} else {
				covStruct@param.n <- as.integer(d)
				covStruct@range.n <- as.integer(d)
			}	
		} else {	
			covStruct@paramset.n <- as.integer(2)
			covStruct@param.n <- as.integer(2*d)
			covStruct@range.n <- as.integer(d)
			covStruct@shape.n <- as.integer(d)
			covStruct@shape.names <- "p"
		}

		if (length(coef.cov)>0) covStruct <- vect2covparam(covStruct, coef.cov)

	} else if (classType=="covAffineScaling") {
   	   	
		covStruct@paramset.n <- as.integer(1)
		covStruct@param.n <- as.integer(2*d)
		covStruct@knots <- c(0,1)
		if (length(coef.cov)>0) covStruct@eta <- coef.cov
    
	} else {
    
    eta.flag <- (length(coef.cov)>0)
		
    for (i in 1:length(knots)) {
      if (is.unsorted(knots[[i]])) {        
        ordKnots <- sort(knots[[i]], index.return = TRUE)
        knots[[i]] <- ordKnots$x
        if (eta.flag) {
          if (length(eta[[i]]) != length(knots[[i]])) stop("mismatch between number of knots and number of values at knots")
          eta[[i]] <- eta[[i]][ordKnots$ix]
        }
      }
    }
    
    names(knots) <- var.names
    covStruct@knots <- knots
    covStruct@param.n <- sum(sapply(knots, length))
    covStruct@paramset.n <- as.integer(1) 
    if (eta.flag) covStruct@eta <- eta
    
	}
  
  validObject(covStruct)
	return(covStruct)	
}

