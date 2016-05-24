# Robust (least absolute deviations LAD/L1) fit of a two-way linear model 
# given a data frame with three columns (x, factor 1, factor 2) or 
# a matrix M where rows of M are ranks of factor 1 levels
# and columns of M are ranks of factor 2 levels, missing values are assumed to be NA or 0.
# No general mean and sum of f2 is constrained to be zero. 
# Estimation method: interior point method in case of dense implementation,
# else Barrodale-Roberts
# A dense matrix implementation is used.
# Some parameters of the dense algorithm are set according to my experience and should be OK
# for most cases. However, they are only tested for sparse design matrices up to the order
# ~ 90.000x2.900
# Calculation is based on treaments contrast for efficiency reason, however,
# they converted to sum contrasts.
# Output: 
#	f1 : parameter estimations of factor 1 (year effects)
#	f1.lev : levels of first factor
#	f2 : parameter estimations of factor 2 (station effects)
#	f2.lev : levels of second factor
#	resid : residuals
#	ierr : return code of the l1 estimation
#	D  : the input as ordered data frame, ordered first after f2 then f1
#	fit  : rq.fit object
pheno.lad.fit <- function(D,limit=1000) {
	if(!is.data.frame(D) && !is.matrix(D)) {
		stop("lad.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(D) && length(D)!=3) {
		stop("lad.fit: argument must be data frame with 3 columns or matrix")
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

	# treatment effects model matrix
	mm <- model.matrix(~ f1 + f2 - 1,na.action=na.exclude)

	if(n > limit) {
		ddm <- as.matrix.csr(mm)
		m <- ddm@dimension[2]
		nnzdmax <- ddm@ia[n + 1] - 1
		l1fit <- rq.fit.sfn(ddm,o,tau=0.5,control=list(tmpmax=1000*m,nnzlmax=100*nnzdmax,small=1e-06))
	}
	else {
		l1fit <- rq.fit(mm,o,tau=0.5,method="br")
	}

	# converting to sum contrasts
	
	# the first station effect
	p2 <-l1fit$coef[-(1:n1)]
	s1 <- -sum(p2)/n2
	p2 <-append(s1,p2+s1)
	# year effects
	p1 <- as.vector(l1fit$coef[1:n1],"numeric")-s1
	resid <- o - (p1[match(f1,levels(f1))]+p2[match(f2,levels(f2))])
	ierr <- l1fit$ierr
	
	return(list(f1=p1,f1.lev=levels(f1),f2=p2,f2.lev=levels(f2),resid=resid,ierr=ierr,D=D,fit=l1fit))
}
