# Fits a two-way linear mixed model given a data frame with three columns
# (x, factor 1, factor 2) or a matrix M where rows of M are ranks of factor 1 levels
# and columns of M are ranks of factor 2 levels, missing values are assumed to be 0.
# The model assumes the first factor f1 to be fixed and the second factor f2 to
# be random. Errors are assumed to be i.i.d. No general mean and sum of f2 is constrained
# to be zero. Estimation method: restricted maximum likelihood (REML)
# Output:
#	fiexed: fixed effects
#	fixed: levels fixed effects
#	random: random effects
#	random.lev: levels random effects
#	resid : residuals
#	SEf1  : standard error group f1 (square root of variance component fixed effect)
#	SEf2  : standard error group f2 (square root of variance component random effect)
#	lclf  : lower 95% confidence limit of fixed effects 
#	uclf  : upper 95% confidence limit of fixed effects 
#	fit  : lme fit object
pheno.mlm.fit <- function(D) {
	if(!is.data.frame(D) && !is.matrix(D)) {
		stop("mlm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.data.frame(D) && length(D)!=3) {
		stop("mlm.fit: argument must be data frame with 3 columns or matrix")
	}
	if(is.matrix(D)) {
		D <- matrix2raw(D)
	}

	D <-  D[order(D[[3]],D[[2]]),]

	s <- factor(D[[3]])
	y <- factor(D[[2]])
	o <- as.vector(D[[1]],"numeric")
	remlfit <- lme(o ~ y - 1 ,random = ~ 1 | s, method="REML", contrasts=list(s=("contr.sum")),na.action=na.exclude)

	fixed <- as.vector(fixed.effects(remlfit),"numeric")
	random <- as.vector(random.effects(remlfit)[[1]],"numeric")
	resid <- as.vector(residuals(remlfit),"numeric")
	SEf1 <-  summary(remlfit)$sigma
	SEf2 <-   attr(remlfit$apVar,"Pars")[[2]]
	lclf <- as.vector(intervals(remlfit,which="fixed")$fixed[,1],"numeric")
	uclf <- as.vector(intervals(remlfit,which="fixed")$fixed[,3],"numeric")
	return(list(fixed=fixed,fixed.lev=levels(y),random=random,random.lev=levels(s),resid=resid,SEf1=SEf1,lclf=lclf,uclf=uclf,D=D,fit=remlfit))
}
