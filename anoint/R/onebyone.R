obo <- function(anoint,...){

	OBO <- onebyone.int(anoint)
	OBO.null <- onebyone.null(anoint)
	LRTs <- LRT.onebyone(OBO,OBO.null)
	names(LRTs) <- all.vars(formula(paste("~",as.character(anoint@formula@prognostic)[3],sep="",collapse="")))
	pvalues <- sapply(LRTs,function(x)pchisq(x,df=1,lower.tail=FALSE))
	names(pvalues) <- names(LRTs)

list(
	fit = OBO,
	LRT = LRTs,
	pvalue = 	pvalues
	)

}

onebyone.int <- function(anoint,...){
	
	if(any(anoint@formula@family=="coxph")){
		lapply(anoint@formula@uni, function(f) coxph(f, data=anoint@data,...))
	}
	else{		
		lapply(anoint@formula@uni, function(f) glm(f, data=anoint@data,
						family=anoint@formula@family,...))
	}
	
}

onebyone.null <- function(anoint,...){
	
	null.formulas <- sapply(anoint@formula@uni,make.null.formula,USE.NAMES=FALSE) 
	
	if(any(anoint@formula@family=="coxph")){
		lapply(null.formulas, function(f) coxph(f, data=anoint@data,...))
	}
	else{		
		lapply(null.formulas, function(f) glm(f, data=anoint@data,
						family=anoint@formula@family,...))
	}
	
}
