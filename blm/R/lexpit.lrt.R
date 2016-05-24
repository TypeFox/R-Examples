LRT <- function(object, var){
	
	if(class(object)[1]!="lexpit"&class(object)[1]!="blm")	
		stop("Object must be of the lexpit or blm class.")
	if(!all(object@weights==1))
		stop("Likelihood-ratio test is only valid for simple random samples.")
	if(class(object)[1]=="lexpit")
		lrt <- LRT.lexpit(object, var)
	else
		lrt <- LRT.blm(object, var)
	
	stat <- 2*(object@loglik-lrt)
	table <- cbind(LRT=stat, pvalue=1-pchisq(stat,df=1))	
	row.names(table) <- var
	
table
}


LRT.blm <- function(object, var){
	
	variable.name.correct <- length(grep(var, attr(terms(object@formula),"term.labels")))!=0
	
	if(!variable.name.correct)
		stop("Variable must be a term label.")
		
	updated.formula <- update(object@formula,paste("~.-",var,sep=""))
	
	fit <- blm(updated.formula, data=object@data, weights=object@weights, strata=object@strata)

fit@loglik
}

LRT.lexpit <- function(object, var){
	
	variable.linear.correct <- length(grep(var, attr(terms(object@formula.linear),"term.labels")))!=0
	variable.expit.correct <- length(grep(var, attr(terms(object@formula.expit),"term.labels")))!=0
	
	if(!variable.linear.correct&!variable.expit.correct)
		stop("Variable must be a term label.")
		
	if(variable.linear.correct){
		formula.linear <- update(object@formula.linear,paste("~.-",var,sep=""))
		formula.expit <- object@formula.expit
	}
	else{
		formula.expit <- update(object@formula.expit,paste("~.-",var,sep=""))
		formula.linear <- object@formula.linear
	}
	
	if(length(attr(terms(formula.linear),"term.labels"))==0){
		fit <- glm(formula.expit, data=object@data, weights=object@weights)
		logLik(fit)
	}
	else if(length(attr(terms(formula.expit),"term.labels"))==0){
		fit <- blm(formula.linear, data=object@data, weights=object@weights)
		fit@loglik
	}
	else{

		fit <- lexpit(formula.linear, formula.expit, data=object@data, weights=object@weights)
					 	
		fit@loglik
	}
}
