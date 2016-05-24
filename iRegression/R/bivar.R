bivar <- function (formula1, lig1, formula2, lig2, data, ...) UseMethod("bivar")

bivar.default <- function (formula1, lig1, formula2, lig2, data, ...) {

	bivarEst <- function (x1, x2, y1, y2, lig1, lig2) {
		
		ro.update = 0
		
		## chutes iniciais
		
		beta1 <- as.vector(coef(lm(y1 ~ x1)))
		beta2 <- as.vector(coef(lm(y2 ~ x2)))
		x1 = cbind(1,x1)
		x2 = cbind(1,x2)
		ro <- as.numeric(cor(y1, y2))
		eta1 <- x1%*%beta1
		eta2 <- x2%*%beta2

		if (lig1 == "identity") 	{mi1 <- eta1; 		deriv_mi1 <- rep(1,length(mi1))} 		else	# identidade
		if (lig1 == "inverse") 		{mi1 <- 1/eta1; 	deriv_mi1 <- as.numeric(-1/(mi1)^2)}		else	# inversa
		if (lig1 == "log") 		{mi1 <- exp(eta1); 	deriv_mi1 <- as.numeric(1/mi1)}			else	# log
		print("Invalid input for the 3th argument. 
		Please enter identity for identity, inverse for inverse and log for logarithmic.")

		if (lig2 == "identity") 	{mi2 <- eta2; 		deriv_mi2 <- rep(1,length(mi2))} 		else	# identidade
		if (lig2 == "inverse") 		{mi2 <- 1/eta2; 	deriv_mi2 <- as.numeric(-1/(mi2)^2)}		else	# inversa
		if (lig2 == "log") 		{mi2 <- exp(eta2); 	deriv_mi2 <- as.numeric(1/mi2)}			else	# log
		print("Invalid input for the 6th argument. 
		Please enter identity for identity, inverse for inverse and log for logarithmic.")

		q <- function(arg1, arg2, ro) {
			return(arg1-(ro*arg2))
		}

		func_c <- function(y1, y2, ro, fi) {
			return((ro*y1*y2 - (y1^2 + y2^2)/2)/fi 
				- log(2*pi*fi/sqrt(1-ro^2))
			)
		}

		b <- function(arg1, arg2, ro) {
			return(
				0.5*((arg1/(1-ro^2) + (arg2*ro)/(1-ro^2))^2 +
				(arg2/(1-ro^2) + (arg1*ro)/(1-ro^2))^2)
				- ro*(arg1/(1-ro^2) + (arg2*ro)/(1-ro^2)) *
				(arg2/(1-ro^2) + (arg1*ro)/(1-ro^2))
			)
		}
		
		di <-	(
			y1 * (q(y1, y2, ro) - q(mi1, mi2, ro)) +
			y2 * (q(y2, y1, ro) - q(mi2, mi1, ro)) +
			b(q(mi1, mi2, ro),q(mi2, mi1, ro),ro) - 
			b(q(y1, y2, ro),q(y2, y1, ro),ro)
		)

		D <- 2*sum(di)
		fi <- D/(2*nrow(x1) - (length(beta1) + length(beta2)))
		ro.seq <- seq(from = -1, to = 1, by=0.001)
		ro.seq <- ro.seq[c(-1,-length(ro.seq))]
		flag = 1
		inter = 0
		par.atual = c(beta1,beta2,ro,fi)

		while (flag != 0) {

			par.atual = c(beta1,beta2,ro,fi)

			# { W

			V1 <- (
				0.5 * (2 * (1/(1 - ro^2) * (1/(1 - ro^2))) + 2 * (ro/(1 - ro^2) *
				(ro/(1 - ro^2)))) - (ro * (1/(1 - ro^2)) * (ro/(1 - ro^2)) +
				ro * (1/(1 - ro^2)) * (ro/(1 - ro^2)))
			)

			V2 <- (	
				0.5 * (2 * (ro/(1 - ro^2) * (ro/(1 - ro^2))) + 2 * (1/(1 - ro^2) *
				(1/(1 - ro^2)))) - (ro * (ro/(1 - ro^2)) * (1/(1 - ro^2)) +
				ro * (ro/(1 - ro^2)) * (1/(1 - ro^2)))
			)

			W1 <- diag(((V1)*(deriv_mi1)^2)^(-1))
			W2 <- diag(((V2)*(deriv_mi2)^2)^(-1))
			zeros <- matrix(0,nrow(x1),nrow(x1))
			W <- cbind(rbind(W1,zeros),rbind(zeros,W2))

			# } W

			# z {

			G1 <- diag(deriv_mi1, nrow = nrow(x1))
			G2 <- diag(deriv_mi2, nrow = nrow(x2))
			z1 <- G1%*%(y1 - mi1)
			z2 <- G2%*%(y2 - mi2)
			z <- rbind(z1,z2)

			# } z

			beta <- rbind(as.matrix(beta1), as.matrix(beta2))

			zerosx <- matrix(0,nrow(x1),ncol(x1))
			X <- cbind(rbind(x1,zerosx),rbind(zerosx,x2))

			betan <- beta + solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%z)

			beta1 <- matrix(betan[seq(from=1,length.out=ncol(x1))])
			beta2 <- matrix(betan[seq(from=1+ncol(x1),length.out=ncol(x1))])

			###############Atualiza Rho ##################

			for (i in 1:length(ro.seq)) {

				log.ver.rho.theta = (y1*q(mi1,mi2,ro.seq[i])+y2*q(mi2,mi1,ro.seq[i]) - 
					b((q(mi1,mi2,ro.seq[i])),(q(mi2,mi1,ro.seq[i])),ro.seq[i]))
				
				log.ver.rho.fi = as.matrix(func_c(y1,y2,ro.seq[i],fi))	

				ro.update[i] = fi^(-1)*sum(log.ver.rho.theta) + sum(log.ver.rho.fi)
	
			}

			ro = ro.seq[which(ro.update[] == max(ro.update))]

			###################################################

			eta1 <- x1%*%beta1
			eta2 <- x2%*%beta2

			if (lig1 == "identity") 	{mi1 <- eta1; 		deriv_mi1 <- rep(1,length(mi1))} 		else	# identidade
			if (lig1 == "inverse") 		{mi1 <- 1/eta1; 	deriv_mi1 <- as.numeric(-1/(mi1)^2)}		else	# inversa
			if (lig1 == "log") 		{mi1 <- exp(eta1); 	deriv_mi1 <- as.numeric(1/mi1)}			else	# log
			print("Invalid input for the 3th argument. 
			Please enter identity for identity, inverse for inverse and log for logarithmic.")

			if (lig2 == "identity") 	{mi2 <- eta2; 		deriv_mi2 <- rep(1,length(mi2))} 		else	# identidade
			if (lig2 == "inverse") 		{mi2 <- 1/eta2; 	deriv_mi2 <- as.numeric(-1/(mi2)^2)}		else	# inversa
			if (lig2 == "log") 		{mi2 <- exp(eta2); 	deriv_mi2 <- as.numeric(1/mi2)}			else	# log
			print("Invalid input for the 6th argument. 
			Please enter identity for identity, inverse for inverse and log for logarithmic.")
	
			di <-	(
				y1 * (q(y1, y2, ro) - q(mi1, mi2, ro)) +
				y2 * (q(y2, y1, ro) - q(mi2, mi1, ro)) +
				b(q(mi1, mi2, ro),q(mi2, mi1, ro),ro) - 
				b(q(y1, y2, ro),q(y2, y1, ro),ro)
			)

			D <- 2*sum(di)

			fi <- D/(2*nrow(x1) - (length(beta1) + length(beta2)))

			par.novo = c(beta1,beta2,ro,fi)
			flag = sum(abs((par.novo-par.atual)/par.atual)>0.0001)
			inter = inter + 1

		} # fim do while

		list(coefficients1 = as.matrix(beta1),
			mi1 = mi1,
			deriv_mi1 = deriv_mi1,
			coefficients2 = as.matrix(beta2),
			mi2 = mi2,
			deriv_mi2 = deriv_mi2,
			inter = inter,
			di = di,
			D = D,
			rho = ro,
			phi = fi)
	}
	
	## extract terms
	mf1 <- model.frame(formula=formula1,data=data)
	x1 <- model.matrix(attr(mf1, "terms"), data=mf1)
	y1 <- model.response(mf1)
	lig1 <- lig1
	mf2 <- model.frame(formula=formula2,data=data)
	x2 <- model.matrix(attr(mf2, "terms"), data=mf2)
	y2 <- model.response(mf2)
	lig2 <- lig2
	## calc
	x1 <- as.matrix(x1[,-1])
	x2 <- as.matrix(x2[,-1])
	y1 <- as.numeric(y1)
	y2 <- as.numeric(y2)
	est <- bivarEst(x1, x2, y1, y2, lig1, lig2)
	
	est$fitted.values1 <- as.vector(est$mi1)
 	est$fitted.values2 <- as.vector(est$mi2)
   	est$residuals1 <- y1 - est$fitted.values1
    	est$residuals2 <- y2 - est$fitted.values2
	est$rd <- as.vector(sign(est$residuals1+est$residuals2)*sqrt(est$di))
   	est$call <- match.call()
    	class(est) <- "bivar"
    	est

}

print.bivar <- function(x, ...) {

	cat("Call:\n")
	print(x$call)
	cat("\n")
	print(list(coefficients1 = x$coefficients1, coefficients2 = x$coefficients2,
		fitted.values1 = x$fitted.values1, fitted.values2 = x$fitted.values2,
		residuals1 = x$residuals1, residuals2 = x$residuals2, residual.deviance = x$rd, 
		rho = x$rho,
		phi = x$phi,
		D = x$D)
	)
}

summary.bivar <- function(object, ...) {

	rmse1 <- sqrt(mean(object$residuals1^2))
	rmse2 <- sqrt(mean(object$residuals2^2))
	res <- list(call = object$call,
		coefficients1 = object$coefficients1,
		RMSE1 = rmse1,
		coefficients2 = object$coefficients2,
		RMSE2 = rmse2,
		rho = object$rho,
		phi = object$phi,
		D = object$D)
	class(res) <- "summary.bivar"
	res
}

print.summary.bivar <- function(x, ...) {

	cat("Call:\n")
	print(x$call)
	cat("\n")
	cat("Coefficients1:\n")
	print(x$coefficients1)
	cat("\n")
	cat("Coefficients2:\n")
	print(x$coefficients2)
	cat("\n")
	cat("RMSE1:\n")
	print(x$RMSE1)
	cat("\n")
	cat("RMSE2:\n")
	print(x$RMSE2)
	cat("\n")
	cat("Rho:\n")
	print(x$rho)
	cat("\n")
	cat("Phi:\n")
	print(x$phi)
	cat("\n")
	cat("D:\n")
	print(x$D)
}

coef.bivar <- function(object, ...) {

	coef1 <- object$coefficients1
	coef2 <- object$coefficients2
	coef <- list(coefficients1 = coef1,
		coefficients2 = coef2)
	class(coef) <- "coef.bivar"
	coef
}

print.coef.bivar <- function(x, ...) {

	print(list(coefficients1 = x$coefficients1,
		coefficients2 = x$coefficients2))
}

fitted.bivar <- function(object, ...) {

	fit1 <- object$fitted.values1
	fit2 <- object$fitted.values2
	ftd <- cbind(fit1,
		fit2)
	fitted <- round(ftd,digits=3)
	class(fitted) <- "fitted.bivar"
	fitted
}

residuals.bivar <- function (object, ...) {

	resid1 <- object$residuals1
	resid2 <- object$residuals2
	resid.dev <- object$rd
	resi <- cbind(resid1,resid2, resid.dev)
	resi = round(resi,digits=3)
	class(resi) <- "residuals.bivar"
	resi
}

bivar.formula <- function(formula1, lig1, formula2, lig2, data=list(),...) {
	
	mf1 <- model.frame(formula=formula1,data=data)
	x1 <- model.matrix(attr(mf1, "terms"), data=mf1)
	y1 <- model.response(mf1)
	mf2 <- model.frame(formula=formula2,data=data)
	x2 <- model.matrix(attr(mf2, "terms"), data=mf2)
	y2 <- model.response(mf2)
	est <- bivar.default(formula1, lig1, formula2, lig2, data, ...)
	est$call <- match.call()
	est$formula1 <- formula1
	est$lig1 <- lig1
	est$formula2 <- formula2
	est$lig2 <- lig2
	est
}