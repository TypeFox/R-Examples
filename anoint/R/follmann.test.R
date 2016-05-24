# FOLLMANN TEST
pim.approx <- function(formula.anoint,data){

follmann.test <- function(control.fit,treated.fit){
		
		coef0 <- coef(control.fit)
        coef1 <- coef(treated.fit)
                
		vmat1 <- vcov(control.fit)
		vmat2 <- vcov(treated.fit)
		vmat <- (vmat1 + vmat2)/2  # Sigma.hat
		# IF GLM, REMOVE INTERCEPT TERM
		if(class(control.fit)[1]!="coxph"){
			coef0 <- coef0[-1]
			coef1 <- coef1[-1]
			vmat <- vmat[2:nrow(vmat),2:nrow(vmat)]	
		}
                
		vc <- solve(t(chol(vmat)))
		
		U0 <- c(vc %*% coef0)
		U1 <- c(vc %*% coef1)
		R <- sqrt(sum(U1*U1)/sum(U0*U0))
		cost <- sum(U0 * U1) / sqrt(sum(U1*U1) * sum(U0 * U0))
		a <- ((R - 1/R) + sqrt((R - 1/R)^2 + 4 * cost^2)) / (2 * cost) # (4)
		u0 <- (a * U1 + U0) / (1 + a^2) # (3)
		u0.norm <- sqrt(sum(u0*u0))/length(u0) 
		u1 <- a * u0 # (5)
		T <- sum((U1 - U0)^2)/2 - sum((U1 - u1)^2) - sum((U0 - u0)^2)

                beta.0 <- solve(vc)%*%u0
        
		list(
				LRT = T,
				theta =a,
				beta.control = beta.0
			)
}

	trt.index <- names(data)==formula.anoint@trt

	if(formula.anoint@family=="coxph"){
		fit0 <- coxph(formula.anoint@prognostic,data[data[,trt.index]==0,])
		fit1 <- coxph(formula.anoint@prognostic,data[data[,trt.index]==1,])			
	}
	else{
		fit0 <- glm(formula.anoint@prognostic,data[data[,trt.index]==0,],
													family=formula.anoint@family)
		fit1 <- glm(formula.anoint@prognostic,data[data[,trt.index]==1,],
													family=formula.anoint@family)			
	}

	fit <- follmann.test(fit0,fit1)
	
	f.trt <- update(formula.anoint@prognostic,paste("~",
		formula.anoint@trt,"offset",sep="+",collapse=""))

	pim.offset <- offset.pim(formula.anoint,data,trt.index,fit$beta.control,fit$theta)
	data$offset <- offset(pim.offset)
	
	if(formula.anoint@family=="coxph"){
		fit$alpha <- coef(coxph(f.trt,data=data))[-2]
	}
	else{
		fit$alpha <- coef(glm(f.trt,data=data,family=formula.anoint@family))[-3]
	}

fit
}
