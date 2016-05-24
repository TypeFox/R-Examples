# Fits a two-way linear model given a data frame with three columns
# (x, factor 1, factor 2) or a matrix M where rows of M are ranks of factor 1 levels
# and columns of M are ranks of factor 2 levels, missing values are assumed to be NA or 0.
# The model assumes the both factors f1 and f2 to be fixed.
# Errors are assumed to be i.i.d. No general mean and sum of f2 is constrained
# to be zero. 
# For n > 1000, for efficiency resaons, the linear model is calcualted
# for treatment contrasts and the constraint that the sum of f2 is zero,
# is adjusted afterwards. This results in a slight over-estimation of
# standard errors.
# Output: f1 : coefficients of first factor
#        f1.lev: levels of first factor
#	 f2 :  coefficients of second fector
#        f2.lev: levels of second factor
#	 resid : residuals
#	 lclf1  : lower 95% confidence limits of factor f1
#	 uclf1  : upper 95% confidence limits of factor f1
#	 lclf2  : lower 95% confidence limits of factor f2
#	 uclf2  : upper 95% confidence limits of factor f2
#	 fit  : the fitted model object
pheno.flm.fit <- function(D,limit=1000) {
	if(!is.data.frame(D) && !is.matrix(D)) {
		stop("flm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(D) && length(D)!=3) {
		stop("flm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.matrix(D)) {
		D <- matrix2raw(D)
	}
	
	
	D <-  D[order(D[[3]],D[[2]]),]

	o <- as.vector(D[[1]],"numeric") # observations
	n <- length(o)	 		# number of observations
	f1 <- factor(D[[2]]) 		# factor 1: year
	n1 <- nlevels(f1) 		# number of levels factor 1 (phenlogy: years)
	f2 <- factor(D[[3]])		# factor 2: station
	n2 <- nlevels(f2)		# number of levels factor 2 (phenlogy: station)

	if(n > limit) { # treatment contrasts
		
		#dense design matrix with treatment contrasts
		ddm <- as.matrix.csr(model.matrix(~ f1 + f2 - 1,na.action=na.exclude))

		m <- ddm@dimension[2]
		fit <- slm.fit(ddm,o,tmpmax=1000*m,small=1e-06)

		fit$terms <- terms(o ~ f1 + f2, na.action=na.exclude)
		fitsum <- summary.slm(fit)

		# converting to sum contrasts
		# the first station effect
		p2 <- as.vector(fitsum$coef[-(1:n1),1],"numeric")
		s1 <- -sum(p2)/n2
		p2 <- append(s1,p2+s1)
		
		p1 <- as.vector(fitsum$coef[1:n1,1],"numeric")-s1	

		# the standard error are not correctly estimated. This has to be adjusted.
		p1.se <- as.vector(fitsum$coef[1:n1,2],"numeric")	

		p2.se <- as.vector(fitsum$coef[-(1:n1),2],"numeric")	
		p2.se <- append(mean(p2.se),p2.se)

		resid <- as.vector(residuals(fit),"numeric")

		df <- fitsum$df[2]

		lclp1 <- p1-qt(0.975,df)*p1.se
		uclp1 <- p1+qt(0.975,df)*p1.se

		lclp2 <- p2-qt(0.975,df)*p2.se
		uclp2 <- p2+qt(0.975,df)*p2.se
	}
	else {		# sum contrasts
		
		#dense design matrix with sum contrasts
		ddm <- pheno.ddm(D)$ddm

		m <- ddm@dimension[2]
		fit <- slm.fit(ddm,append(o,0),tmpmax=1000*m,small=1e-06)

		fit$terms <- terms(o ~ f1 + f2, na.action=na.exclude)
		fitsum <- summary.slm(fit)

		p1 <- as.vector(fitsum$coef[1:n1,1],"numeric")	
		p1.se <- as.vector(fitsum$coef[1:n1,2],"numeric")	

		p2 <- as.vector(fitsum$coef[-(1:n1),1],"numeric")	
		p2.se <- as.vector(fitsum$coef[-(1:n1),2],"numeric")	

		resid <- as.vector(residuals(fit),"numeric")

		df <- fitsum$df[2]

		lclp1 <- p1-qt(0.975,df)*p1.se
		uclp1 <- p1+qt(0.975,df)*p1.se

		lclp2 <- p2-qt(0.975,df)*p2.se
		uclp2 <- p2+qt(0.975,df)*p2.se
	}

	return(list(f1=p1,f1.se=p1.se,f1.lev=levels(f1),f2=p2,f2.se=p2.se,f2.lev=levels(f2),resid=resid,lclf1=lclp1,uclf1=uclp1,lclf2=lclp2,uclf2=uclp2,D=D,fit=fit))
}
