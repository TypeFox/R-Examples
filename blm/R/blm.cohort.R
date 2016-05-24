blm.cohort <- function(formula,data,na.action=na.omit,par.init=NULL,warn=FALSE,...){

#### DEPENDENT FUNCTIONS
logit <- function(x) log(x/(1-x))

na.lexpit <- function(f,data,FUN){ # REMOVING MISSING
	
	keep.data <- subset(data,select=all.vars(f))
	keep.data <- FUN(keep.data)
	kept <- match(row.names(keep.data),row.names(data))

list(data = keep.data, kept = kept)
}

LL <- function(Y,p){ # LOG-LIKELIHOOD
 -sum(Y*logit(p)+log(1-p))
}

blm.cohort.optim <- function(Y,X,beta.init,...){ # CONSTRAINED OPTMIZER

loglik <- function(beta,Y,X){
	
	p <- X%*%beta
	
LL(Y,p)
}

gradient <- function(beta,Y,X){
	
	p <- X%*%beta
	C <- (Y*1/(p*(1-p))-1/(1-p))
	C <- matrix(C,nrow(X),ncol(X))

-colSums(C*X)
}

# DEFINE CONSTRAINTS
lexpit.constraints <- function(X){
		
	# CONSTRAINTS OF FORM
	# X%*%beta-c >= 0
	list(
		U = rbind(X,-X),
		C = c(rep(0,nrow(X)),rep(-1,nrow(X)))	
	)
}

	constr <- lexpit.constraints(X)

	results <- constrOptim(beta.init,f=loglik,grad=gradient,
						Y = Y, X = X,ui = constr$U,ci = constr$C,...)

results
}

#####

	warn.option <- getOption("warn")
	if(!warn) options("warn"=-1)
	
	data <- na.lexpit(formula,data,FUN=na.action)
	which.kept <- data$kept
	data <- data$data
	
	Y <- model.frame(formula,data)[,1]
	X <- model.matrix(formula,data)
			
	if(is.null(par.init)){
		beta.init <- rep(0,ncol(X))
		beta.init[1] <- mean(Y)
	}
	else{
		beta.init <- par.init
		}
		
	fit <- blm.cohort.optim(Y,X,beta.init,...)
	beta <- fit$par
				
	names(beta) <- colnames(X)
	
	# NULL MODEL, ESTIMATE IS THE OVERALL MEAN
	ll.null <- -LL(Y,mean(Y))
	
	options("warn"=warn.option)
	
	list(
		coef = beta,
		formula= formula,
		df.residual = nrow(X)-ncol(X),
		data = data,
		converged = fit$convergence==0,
		par.init = beta.init,
		loglik = -LL(Y,X%*%beta),
		loglik.null = ll.null,
		barrier.value = fit$barrier.value
	)
	
}


LRT <- function(object,var=NULL,...){
	
	f <- object$formula
	labels <- attr(terms(f),"term.labels")
	
	if(!is.null(var)) labels <- labels[grep(var,labels)]
	
	formulas <- lapply(labels, function(term) update(f, paste("~.-",term, sep="")))
	
	LL <- sapply(formulas, function(new.formula) blm.cohort(new.formula,object$data)$loglik)
	
	LRTs <- 2*(object$loglik-LL)
	
	table <- cbind(LRT=LRTs, pvalue=1-pchisq(LRTs,df=1))
	row.names(table) <- labels

table
}


