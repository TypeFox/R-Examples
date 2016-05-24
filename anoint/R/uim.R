uim <- function(object,...){
	
	fit.int <- uim.int(object)
	fit.null <- uim.null(object)
	LRT <- LRT.uim(fit.int,fit.null)	
	
list(
	fit = fit.int,
	LRT = LRT,
	pvalue = pchisq(LRT,df=length(object@formula@uni),lower.tail=FALSE)
)
}

uim.int <- function(object,...){
			
	if(object@formula@family=="coxph"){
		uim <- coxph(object@formula@formula,object@data,...)
	}
	else{
		uim <- glm(object@formula@formula,object@data,family=object@formula@family,...)
	}

uim
}

uim.null <- function(object,...){
			
	f <- make.null.formula(object@formula@formula)
	
	if(object@formula@family=="coxph"){
		uim <- coxph(f,object@data,...)
	}
	else{
		uim <- glm(f,object@data,family=object@formula@family,...)
	}

uim
}
