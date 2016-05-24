blm <- function(formula,data,na.action=na.omit,
							weights=NULL,strata=NULL,par.init=NULL,warn=FALSE,...){

na.lexpit <- function(f,data,FUN){

	keep.data <- subset(data,select=all.vars(f))
	keep.data <- FUN(keep.data)
	kept <- match(row.names(keep.data),row.names(data))

list(data = keep.data, kept = kept)
}

LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}

	warn.option <- getOption("warn")
	if(!warn) options("warn"=-1)
	
	data <- na.lexpit(formula,data,FUN=na.action)
	which.kept <- data$kept
	data <- data$data
	
	Y <- model.frame(formula,data)[,1]
	X <- model.matrix(formula,data)
	
	if(is.matrix(X)){
		x.labels <- colnames(X)
	}
	else{
		x.labels <- attr(terms(formula), "term.labels") 
	}
	
	if(is.null(strata)){
				strata <- rep(1,nrow(data)) 
			}
	else{
			strata <- strata[which.kept]
		}
			
	if(!class(strata)=="factor") strata <- factor(strata)
	
	if(is.null(weights)){
			weights <- rep(1,nrow(X))
			}
	else{
			weights <- weights[which.kept]
			}
	
	# STANDARDIZE WEIGHTS FOR STABILITY
	w <- cbind(weights/mean(weights))
			
	if(is.null(par.init)){
		beta.init <- rep(0,ncol(X))
		beta.init[1] <- sum(Y*w)/sum(w)
	}
	else{
		beta.init <- par.init
		}
		
	fit <- blm.optim(Y,X,w,beta.init,...)
	beta <- fit$par
	
	# if(!all(weights==1))
		# vcov <- solve(vcov.blm.revised.strata(beta,cbind(Y),X,cbind(w),strata))
	# else
		# vcov <- solve(vcov.blm.revised.bigmatrix(beta,cbind(Y),X,cbind(w)))
	
	 if(!all(weights==1)){
	 	if(nrow(data)>50000)
		 	vcov <- vcov.blm.big(formula, data, weights, beta)
		 else
		 	vcov <- vcov.influence.blm.strata(formula, data, weights, strata, beta)
		 	}
	 else{
	 	if(nrow(data)>50000)
		 	vcov <- vcov.blm.big(formula, data, weights=NULL, beta)
		 else
		 	vcov <- vcov.influence.blm(formula, data, beta)
	  }

	 	
	names(beta) <- x.labels
	
	# NULL MODEL, ESTIMATE IS THE OVERALL MEAN
	ll.null <- -LL(Y,sum(Y*w)/sum(w),cbind(weights))
	
	options("warn"=warn.option)
	
	new("blm",
		coef = beta,
		vcov = vcov,
		formula= formula,
		df.residual = nrow(X)-ncol(X),
		data = data,
		which.kept = which.kept,
		y = Y,
		weights = as.numeric(weights),
		strata = strata,
		converged = fit$convergence==0,
		par.init = beta.init,
		loglik = -LL(Y,X%*%beta,cbind(weights)),
		loglik.null = ll.null,
		barrier.value = fit$barrier.value
	)
	
}

