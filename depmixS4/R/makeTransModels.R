makeTransModels <-
function(nstates,formula=~1,data=NULL,homogeneous,values=NULL, ...) {
	
	# defaults that possibly need some work at some point 
	# FIX ME
	base=1
	prob=TRUE
	
	if(!homogeneous&is.null(data)) stop("non-homogeneous transition models needs data argument")
	
	# starting values	
	tst <- FALSE
	if(!is.null(values)) {
		tst <- TRUE
		values <- matrix(values,nstates,byrow=TRUE)
	}
	
	models <- list()
	for(i in 1:nstates) {
		if(tst) {
			if(homogeneous) models[[i]] <- transInit(formula,multinomial(link="identity"),data=data.frame(1),nstates=nstates,pstart=values[i,],prob=prob)
			else models[[i]] <- transInit(formula,multinomial(base=base),data=data,nstates=nstates,pstart=values[i,],prob=prob)
		} else {
			if(homogeneous) models[[i]] <- transInit(formula,multinomial(link="identity"),data=data.frame(1),nstates=nstates,prob=FALSE)
			else models[[i]] <- transInit(formula,multinomial(base=base),data=data,nstates=nstates,prob=FALSE)
		}
	}
	
	return(models)
}

