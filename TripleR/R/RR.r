


#' @export
RR.summary <- function(formule, data) {
	# save warnings for later output
	l <- long2matrix(formule, data)
	cat(paste("Number of groups:",length(l))); cat("\n")
	cat(paste("Number of valid groups (i.e., more than 3 group members):\n",sum(laply(l, function(x) ncol(x) > 3)))); cat("\n")
	
	# which groups are excluded?
	excl.id <- c()
	for (i in 1: length(l)) {if (nrow(l[[i]]) <=3) excl.id <- c(excl.id, attr(l[[i]], "group.id"))}
	if (length(excl.id>0)) {
		cat("Following groups are excluded because they have 3 or less members:")
		cat(paste("\n\t",excl.id)) 
		cat("\n\n")
	}

	cat("Group sizes:\n")
	for (i in 1:length(l)) {cat(paste(attr(l[[i]], "group.id"),": n=",ncol(l[[i]]),"\n", sep=""))}

	cat("\n\n")
	cat("Group member statistics:"); cat("\n")
	gs <- laply(l, function(x) ncol(x))	# group sizes
	gm <- table(gs)
	names(gm) <- paste("n=",names(gm), sep="")
	print(gm); cat("\n")
	cat(paste("Min:",min(gs))); cat("\n")
	cat(paste("Max:",max(gs))); cat("\n")
	cat(paste("Mean:",round(mean(gs), 2))); cat("\n\n")
}




# berechnet schnell die Effekte, ohne sonstigen Krimskrams
quickeffects <- function(RRMatrix) {
	RRMatrix <- as.matrix(RRMatrix)
	mip <- rowMeans(RRMatrix, na.rm=TRUE)
	mpj <- colMeans(RRMatrix, na.rm=TRUE)
	mpp <- mean(RRMatrix, na.rm=TRUE)  # grand mean

   n <- length(mip)
   a <- ((n-1)^2)/(n*(n-2))*mip + (n-1)/(n*(n-2))*mpj - (n-1)/(n-2)*mpp #actor effects
   b <- ((n-1)^2)/(n*(n-2))*mpj + (n-1)/(n*(n-2))*mip - (n-1)/(n-2)*mpp #partner effects
   
   am <- matrix(a, nrow(RRMatrix), ncol(RRMatrix))
   bm <- t(matrix(b, nrow(RRMatrix), ncol(RRMatrix)))
   c <- RRMatrix - am - bm - mpp		 # relationship effect

	return(list(a=a,b=b,c=c,m=mpp))
}




# calculates Actor-, Partner- and Relationship-Effects from a single RR-Matrix
RR.effects <- function(RRMatrix, name=NA, na.rm=FALSE, index="", varname="NA") {
	if (!is.na(varname)) {name <- varname} else {
		if (!is.null(attr(RRMatrix, "varname"))) name <- attr(RRMatrix, "varname")
	}
	
	RRold <- RRMatrix		
	if (na.rm & sum(is.na(RRMatrix))>nrow(RRMatrix)) {
		imp <- impute(RRMatrix)
		RRMatrix <- imp$RRMatrix
		imputation <- TRUE
	} else {
		imputation <- FALSE
	}
	

	RRMatrix2 <- as.matrix(RRMatrix)
	mip <- rowMeans(RRMatrix2, na.rm=TRUE)
	mpj <- colMeans(RRMatrix2, na.rm=TRUE)
	mpp <- mean(RRMatrix2, na.rm=TRUE)  # grand mean

	n <- length(mip)
	a <- ((n-1)^2)/(n*(n-2))*mip + (n-1)/(n*(n-2))*mpj - (n-1)/(n-2)*mpp #actor effects
	b <- ((n-1)^2)/(n*(n-2))*mpj + (n-1)/(n*(n-2))*mip - (n-1)/(n-2)*mpp #partner effects

	am <- matrix(a, nrow(RRMatrix2), ncol(RRMatrix2))
	bm <- t(matrix(b, nrow(RRMatrix2), ncol(RRMatrix2)))
	c <- RRMatrix2 - am - bm - mpp		 # relationship effect
	rownames(c) <- colnames(c) <- rownames(RRMatrix)
	
	# delete all relationship effects which had a NA in the original matrix
	if (imputation) {
		c[imp$NAs==TRUE] <- NA
	}
	


	# actor and partner effects
	self <- FALSE
	
	if (!is.null(attr(RRMatrix, "self.ratings"))) {
		self <- TRUE
		self.centered <- attr(RRMatrix, "self.ratings")-mean(attr(RRMatrix, "self.ratings"), na.rm=TRUE)
		
		eff <- data.frame(id = rownames(RRMatrix), actor=a, partner=b, self=self.centered)
		if (!is.null(name)) {colnames(eff)[2:4] <- paste(name, localOptions$suffixes, sep="")}
		attr(eff[,4], "type") <- "self"
		
		## Add selfenhancement index?
		if (index != "") {
			if (match.arg(index, c("enhance"))=="enhance") {
				kwan <- self.centered - a - b
				eff <- data.frame(eff, enhance=kwan)
				if (!is.null(name)) {colnames(eff)[colnames(eff)=="enhance"] <- paste(name, ".enhance", sep="")}
			}
		}
		
	} else {
		eff <- data.frame(id = rownames(RRMatrix), actor=a, partner=b)
		if (!is.null(name)) {colnames(eff)[2:3] <- paste(name, localOptions$suffixes[1:2], sep="")}
	}
	
	attr(eff[,2], "type") <- "actor"
	attr(eff[,3], "type") <- "partner"
	
	
	
	## Relationship effects
	
	effRel <- reshape2::melt(c)
	effRel <- effRel[apply(effRel, 1, function(x) {x[1] != x[2]}),]
	colnames(effRel) <- c("actor.id", "partner.id", "relationship")
	effRel[,1] <- factor(effRel[,1])
	effRel[,2] <- factor(effRel[,2])
	
	# sort relEffects according to dyad
	digits <- floor(log10(n))+1
	effRel$dyad <- factor(apply(effRel, 1, function(x) paste(sort(x[1:2], decreasing=FALSE), collapse="_")))
	effRel$dyad <- factor(effRel$dyad, labels=paste(attr(RRold, "group.id"), "_", f2(1:length(levels(effRel$dyad)), 0, paste("0",digits,sep="")), sep=""))
	effRel <- effRel[,c(1,2,4,3)]
	effRel <- effRel[order(effRel$dyad, effRel$actor.id),]
	
	
	
	## group means added
	
	eff.gm <- eff
	if (!is.null(attr(RRMatrix, "self.ratings"))) {
		eff.gm[,2:3] <- eff.gm[,2:3]+mpp
		eff.gm[,4] <- attr(RRMatrix, "self.ratings")
	} else {
		eff.gm[,2:3] <- eff.gm[,2:3]+mpp
	}
	


	## construct return object
	if (!is.null(attr(RRMatrix, "self.ratings"))) {
		res <- list(actor = a, partner = b, relationship = c, eff=eff, effRel=effRel, eff.gm=eff.gm, self=self.centered)
	} else {
		res <- list(actor = a, partner = b, relationship = c, eff=eff, effRel=effRel, eff.gm=eff.gm)
	}
	
	attr(res, "self") <- self
	
	return(res)
}




# calculates variance components from a single RR-Matrix
RR.univariate <- function(RRMatrix, na.rm=FALSE, verbose=TRUE, corr.fac="1", index="", varname=NA) {
	
	if (is.null(RRMatrix)) return();
	
	if (nrow(RRMatrix)<4) {
		warning(paste("WARNING: group",attr(RRMatrix, "group.id"),"has 3 or fewer subjects. For calculation of SRM variables minimum group size is 4."), call.=FALSE);
		return();
	}
	
	
	# emit some warnings about missings if there are NAs outside the diagonale
	if ((sum(is.na(RRMatrix)) > nrow(RRMatrix)) & na.rm==FALSE)
		stop("There are NAs outside the diagonale. Calculations are aborted.")
		
	n <- nrow(RRMatrix)	
	n.NA <- sum(is.na(RRMatrix)) - n
	
	
	## Warning if too many missings are present
	if (n.NA > 0 & na.rm==TRUE & verbose==TRUE) {
		wa <- FALSE
		if (n==4 & n.NA > 1) {wa <- TRUE}
		if (n==5 & n.NA > 2) {wa <- TRUE}
		if (n==6 & n.NA > 4) {wa <- TRUE}
		if (n==7 & n.NA > 6) {wa <- TRUE}
		if (n==8 & n.NA > 8) {wa <- TRUE}
		if (n>=10 & (n.NA/(n^2-n)) > .20) {wa <- TRUE}
		
		if (wa==TRUE) {
			warning(paste(attr(RRMatrix, "varname"),": The number of missing values (n.NA=",n.NA,"; ",round((n.NA/(n^2-n))*100, 1),"%) in group ",attr(RRMatrix, "group.id")," exceeds the recommended maximum number of missings according to Schoenbrodt, Back, & Schmukle (2011). Estimates might be biased.", sep=""), call.=FALSE)
		}
	}
	
	eff <- RR.effects(RRMatrix, name=attr(RRMatrix, "varname"), na.rm=na.rm, index=index, varname=varname)
	
	A <- sum(eff$actor^2)/(n-1)
	B <- sum(eff$partner^2)/(n-1)
	C <- sum(eff$actor*eff$partner)/(n-1)
	e <- 0.5*(eff$relationship + t(eff$relationship))
	d <- eff$relationship - t(eff$relationship)
	
	if (na.rm==TRUE) {
		
		if (is.na(corr.fac)) {
			corr.fac <- (n*(n-1)) / (n*(n-1) - sum(is.na(e)) + n)
		} else {
			corr.fac <- eval(parse(text=corr.fac))
		}
		
		D <- (sum(e^2, na.rm=TRUE) * corr.fac)  / (((n-1)*(n-2)/2)-1)
		E <- ((sum(d^2,na.rm=TRUE)/2)  * corr.fac) / ((n-1)*(n-2))
	} else {
		D <- sum(e^2, na.rm=TRUE) / (((n-1)*(n-2)/2)-1)
		E <- (sum(d^2,na.rm=TRUE)/2) / ((n-1)*(n-2))
	}
	
	
	scc <- (D+E)/2 #relationship variance
	sccs <- (D-E)/2 #relationship covariance
	sab <- C - (sccs*(n-1))/(n*(n-2)) - scc/(n*(n-2)) #actor-partner covariance
	saa <- A - (scc*(n-1))/(n*(n-2)) - sccs/(n*(n-2)) #actor variance
	sbb <- B - (scc*(n-1))/(n*(n-2)) - sccs/(n*(n-2)) #partner variance
	
	saa2 <- ifelse(saa>=0, saa, NaN)
	sbb2 <- ifelse(sbb>=0, sbb, NaN)
	scc2 <- ifelse(scc>=0, scc, NaN)
	
	raa <- saa2/sum(saa2,sbb2,scc2,na.rm=TRUE) #standardized actor variance
	rbb <- sbb2/sum(saa2,sbb2,scc2,na.rm=TRUE) #standardized partner variance
	rcc <- scc2/sum(saa2,sbb2,scc2,na.rm=TRUE) #standardized relationship variance
	rab <- ifelse(saa>0 & sbb>0,sab/sqrt(saa*sbb),NaN) #actor-partner correlation
	rccs <- sccs/scc2 #relationship correlation
	
	
	# ---------------------------------------------------------------------
	# Compute SE and t-value for a single group using the formula of Lashley & Bond
	SEVAR <- compute_univariate_LB_SE2(saa, sbb, scc, sab, sccs, n) # squared standard errors of variance estimates
	SE <- sqrt(SEVAR)	# standard errors of variance estimates
	
	# error variance is NA if only one group is present
	estimate <- c(saa,sbb,scc,NA,sab,sccs)
	standardized <- clamp(raa,rbb,rcc,NA,rab,rccs)
	
	# set SEs of negative variances to NaN
	# TODO: If var==0 --> se = 0, too?
	SE[estimate[1:3]<=0] <- NaN
	
	t.value <- estimate/SE
	p.value <- 1-pt(abs(t.value), n-1)
	# Kovarianzen werden zweiseitig getestet:
	p.value[4:5] <- p.value[4:5]*2
	
	
	# calculate reliability for actor and partner effects
	rel.a <- saa / (saa + scc*(n-1)/(n*(n-2)) + sccs/(n*(n-2)))
	if (saa < 0) rel.a <- NaN
	rel.b <- sbb / (sbb + scc*(n-1)/(n*(n-2)) + sccs/(n*(n-2)))
	if (sbb < 0) rel.b <- NaN
	
	attr(eff$eff[,2], "reliability") <- rel.a
	attr(eff$eff[,3], "reliability") <- rel.b
	
	# join everything in one dataframe
	univariate <- data.frame(type=unilabels_b, estimate, standardized, se=SE, SEVAR=SEVAR, t.value, p.value)
	#univariate <- data.frame(type=unilabels_b, estimate, standardized, se=SE, t.value, p.value)
	
	# if one variance component is below zero: erase covariances
	# erase indices for negative variances
	univariate[1:3,][univariate$estimate[1:3]<0, c("standardized", "se", "t.value", "p.value")] <- NaN
	if (saa <= 0 | sbb <= 0) {univariate[5, c("standardized", "se", "t.value", "p.value")] <- NaN}
	if (scc <= 0) {univariate[6, c("standardized", "se", "t.value", "p.value")] <- NaN}

	res <- list(effects = eff$eff, effectsRel = eff$effRel, effects.gm = eff$eff.gm, varComp = univariate, relMat.av=e, relMat.diff=d, group.size=n, latent=FALSE, anal.type="Univariate analysis of one round robin variable", n.NA = n.NA, SEVAR=SEVAR)
	class(res) <- "RRuni"
	attr(res, "group.size") <- n
	attr(res, "varname") <- attr(RRMatrix, "varname")
	attr(res, "self") <- attr(eff, "self")
	
	# if self ratings are present: add to results object
	#print(attr(RRMatrix, "group.id"))
	dummy <- capture.output(self <- invisible(selfCor(res)))
	if (!is.null(self)) {res[["selfCor"]] <- self}
	
	return(res)
}



# combines two RR-matrices, depending on parameter 'latent'
# latent = TRUE: both matrices are treated as two measures for one underlying construct
# latent = FALSE: both matrices are treated as independent variables
# noCorrection = TRUE: even if univariate estimates are negative, bivariate covariances are NOT set to NA (this is necessary, when the manifest bivariat results are transferred into the bivariate latent analysis, see TAG1)

RR.bivariate <- function(RRMatrix1, RRMatrix2, analysis="manifest", na.rm=FALSE, verbose=TRUE, noCorrection=FALSE, index="", varname=NA, se="LashleyBond") {
	
	if (!(analysis %in% c("latent", "manifest"))) stop("Parameter 'analysis' must either be 'latent' or 'manifest'. Calculations aborted.")
	
	
	dif1 <- union(setdiff(rownames(RRMatrix1), rownames(RRMatrix2)), setdiff(rownames(RRMatrix2), rownames(RRMatrix1)))
	if (length(dif1)>0 & verbose==TRUE) {
		warning(paste(length(dif1),"participant(s) have been excluded from the bivariate/latent analysis due to missing data in one of both variables"), call.=FALSE)
	}
	
	
	# Clean up data: only participants are allowed that are in BOTH matrices
	allparticipants <- intersect(rownames(RRMatrix1), rownames(RRMatrix2))
	
	a1 <- attributes(RRMatrix1)
	a1$self.ratings <- a1$self.ratings[allparticipants]
	a1$dimnames <- NULL
	a2 <- attributes(RRMatrix2)
	a2$self.ratings <- a2$self.ratings[allparticipants]
	a2$dimnames <- NULL
	a1$dim <- rep(length(allparticipants), 2)
	a2$dim <- rep(length(allparticipants), 2)

	RRMatrix1 <- RRMatrix1[allparticipants,allparticipants]
	RRMatrix2 <- RRMatrix2[allparticipants,allparticipants]
	dimn <- rownames(RRMatrix1)
	attributes(RRMatrix1) <- a1
	attributes(RRMatrix2) <- a2
	rownames(RRMatrix1) <- rownames(RRMatrix2) <- colnames(RRMatrix1) <- colnames(RRMatrix2) <- dimn
	
	
	RR.1 <- RR.univariate(RRMatrix1, na.rm, verbose, index=index, varname=varname)
	RR.2 <- RR.univariate(RRMatrix2, na.rm, verbose, index=index, varname=varname)	
	varComp.1 <- RR.1$varComp$estimate
	varComp.2 <- RR.2$varComp$estimate
	n <- nrow(RRMatrix1)
	
	

	#Bivariate Relationships
	A <- sum(RR.1$effects[,2]*RR.2$effects[,2])/(n-1)
	B <- sum(RR.1$effects[,2]*RR.2$effects[,3])/(n-1)
	C <- sum(RR.1$effects[,3]*RR.2$effects[,2])/(n-1)
	D <- sum(RR.1$effects[,3]*RR.2$effects[,3])/(n-1)
	E <- sum(RR.1$relMat.av*RR.2$relMat.av,na.rm=TRUE)/(((n-1)*(n-2)/2)-1)
	F <- (sum(RR.1$relMat.diff*RR.2$relMat.diff,na.rm=TRUE)/2)/((n-1)*(n-2))
	sch <- (E+F)/2 #intrapersonal relationship covariance
	schs <- (E-F)/2 #interpersonal relationship covariance
	saf <- A - (sch*(n-1))/(n*(n-2)) - schs/(n*(n-2)) #actor-actor covariance
	sag <- B - (schs*(n-1))/(n*(n-2)) - sch/(n*(n-2)) #actor-partner covariance
	sbf <- C - (schs*(n-1))/(n*(n-2)) - sch/(n*(n-2)) #partner-actor covariance
	sbg <- D - (sch*(n-1))/(n*(n-2)) - schs/(n*(n-2)) #partner-partner covariance
	
	
	# standardized covariances (=correlations), bivariate case
	# standardized <- clamp(raf,rbg,rag,rbf,rch,rchs)
	w <- getOption("warn")
	options(warn=-1)
		raf <-  saf/(sqrt(varComp.1[1])*sqrt(varComp.2[1])) # bivariate correlations
		rbg <-  sbg/(sqrt(varComp.1[2])*sqrt(varComp.2[2]))
		rag <-  sag/(sqrt(varComp.1[1])*sqrt(varComp.2[2]))
		rbf <-  sbf/(sqrt(varComp.1[2])*sqrt(varComp.2[1]))
		rch <-  sch/(sqrt(varComp.1[3])*sqrt(varComp.2[3]))
		rchs <- schs/(sqrt(varComp.1[3])*sqrt(varComp.2[3]))
	options(warn=w)
	
	
	if (analysis=="latent") {
		stabpervar1 <- saf
		stabtarvar1 <- sbg
		stabrelvar1 <- sch
		
    	stabapcov1 <- (sag + sbf)/2	# latent actor-partner-covariance
    	stabdycov1 <- schs
		unstabper1 <- (varComp.1[1] + varComp.2[1])/2 - saf
		unstabtar1 <- (varComp.1[2] + varComp.2[2])/2 - sbg
		unstabrel1 <- (varComp.1[3] + varComp.2[3]) / 2 - sch		
		
		saf2 <- max(saf, 0)
		sbg2 <- max(sbg, 0)
		sch2 <- max(sch, 0)
		
		stable1 <- saf2 + sbg2 + sch2
		unstable1 <- max(unstabper1, 0) + max(unstabtar1, 0) + max (unstabrel1, 0) 
		unstable.raw <- sum(unstabper1, unstabtar1, unstabrel1)
		stabler1 <- stable1 / (stable1 + unstable1)
		unstabler1 <- unstable1 / (stable1 + unstable1)
		stabperr1 <- saf2/(stable1 + unstable1)
		stabtarr1 <- sbg2/(stable1+unstable1)
		stabrelr1 <- sch2/(stable1+unstable1)
		stabdycor1 <- ifelse(sch2>0, stabdycov1/sch2, NaN)
		stabapcor1 <- ifelse(saf>0 & sbg>0, stabapcov1 / sqrt(saf*sbg), NaN)
	}

	# Compute standard errors (se) und t-values (t) of bivariate srm-parameters
	biSEVAR <- compute_bivariate_LB_SE2(varComp.1, varComp.2, saf, sag, sbg, sbf, sch, schs, n)
	
	biSE <- sqrt(biSEVAR)
	names(biSE) <- c("sesaf", "sesbg", "sesag", "sesbf", "sesch", "seschs")

	taf <- saf/biSE["sesaf"]
	tbg <- sbg/biSE["sesbg"]
	tch <- sch/biSE["sesch"]
	tag <- sag/biSE["sesag"]
	tbf <- sbf/biSE["sesbf"]
	tchs <- schs/biSE["seschs"]
	
	if (analysis=="latent") {
		sestabpervar1 <- biSE["sesaf"]
		sestabtarvar1 <- biSE["sesbg"]
		sestabrelvar1 <- biSE["sesch"]

		tstabpervar1 <- ifelse(saf>=0, saf/biSE["sesaf"], NaN)
		tstabtarvar1 <- ifelse(sbg>=0, sbg/biSE["sesbg"], NaN)
		tstabrelvar1 <- ifelse(sch>=0, sch/biSE["sesch"], NaN)
	}

	#########################################Result Matrix


# Result Matrix for two independent variables

if (analysis=="manifest") {

	univariate <- list()
	univariate[[1]] <- RR.1
	univariate[[2]] <- RR.2
	
	estimate <- c(saf,sbg,sag,sbf,sch,schs)
	standardized <- clamp(raf,rbg,rag,rbf,rch,rchs)
	
	t.value <- c(taf, tbg, tag, tbf, tch, tchs)
	pvalues <- (1-pt(abs(t.value), n-1))*2 	# alles Kovarianzen, daher zweiseitig testen!
	bivariate <- data.frame(type=bilabels_bb, estimate, standardized, se=biSE, biSEVAR=biSEVAR, t.value, p.value=pvalues)
	#bivariate <- data.frame(type=bilabels_bb, estimate, standardized, se=biSE, t.value, p.value=pvalues)
	
	if (noCorrection==FALSE) {
		# erase covariances if one variance component is < 0
		if (RR.1$varComp[1,2] <= 0) bivariate[c(1,3), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (RR.1$varComp[2,2] <= 0) bivariate[c(2,4), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (RR.2$varComp[1,2] <= 0) bivariate[c(1,4), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (RR.2$varComp[2,2] <= 0) bivariate[c(2,3), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (RR.1$varComp[3,2] <= 0) bivariate[c(5,6), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (RR.2$varComp[3,2] <= 0) bivariate[c(5,6), c("standardized", "se", "t.value", "p.value")] <- NaN
	}	
	
	res <- list(univariate = univariate, bivariate = bivariate, latent=FALSE, anal.type="Bivariate analysis of two variables, each measured by one round robin variable", biSEVAR=biSEVAR)
	attr(res, "group.size") <- n
	class(res) <- "RRbi"
	
} else 
{
	# Result matrix for latent analysis 

	unstand <- c(stabpervar1,stabtarvar1,stabrelvar1,unstable1,stabapcov1,stabdycov1)
	stand <- clamp(stabperr1, stabtarr1, stabrelr1, unstabler1,stabapcor1,stabdycor1)
	stand[is.infinite(stand)] <- NaN
	
	# se, t, und p werden aus dem bivariaten Fall uebernommen
	se <- c(sestabpervar1, sestabtarvar1, sestabrelvar1, NA, sqrt((biSE["sesag"]^2+biSE["sesbf"]^2)/2), biSE["seschs"])
	tvalues <- c(tstabpervar1,tstabtarvar1,tstabrelvar1,NA,stabapcov1/sqrt((biSE["sesag"]^2+biSE["sesbf"]^2)/2), tchs)
	pvalues <- (1-pt(abs(tvalues), n-1))
	pvalues[4:5] <- pvalues[4:5]*2
	
	results <- data.frame(type=unilabels_b, estimate=unstand, standardized=stand, se=se, SEVAR=c(biSEVAR[1:2], biSEVAR["sesch2"], NA, (biSE["sesag"]^2+biSE["sesbf"]^2)/2, biSEVAR[6]), t.value=tvalues, p.value=pvalues)
	#results <- data.frame(type=unilabels_b, estimate=unstand, standardized=stand, se=se, t.value=tvalues, p.value=pvalues)
	
	# erase indices for negative variances
	results[1:3,][results$estimate[1:3]<0, c("standardized", "se", "t.value", "p.value")] <- NaN
	
	#-------------------------------
	# calculate reliability for actor and partner effects
	
	unstabdycov1 <- (varComp.1[5] + varComp.2[5]) / 2 - schs
	
	r <- 2 # r = number of replications - in our case, it's always 2
	rel.a <- stabpervar1 / ((stabpervar1 + (unstabper1/r)) + (stabrelvar1+(unstabrel1/r))*(n-1)/(n*(n-2)) + (stabdycov1+(unstabdycov1/r))/(n*(n-2)))
	if (stabpervar1 < 0) rel.a <- NaN
	
	rel.p <- stabtarvar1 / ((stabtarvar1 + (unstabtar1/r)) + (stabrelvar1+(unstabrel1/r))*(n-1)/(n*(n-2)) + (stabdycov1+(unstabdycov1/r))/(n*(n-2)))
	if (stabtarvar1 < 0) rel.p <- NaN
	
	rel.r <- stabrelvar1 / (stabrelvar1+(unstabrel1/r))
	if (stabrelvar1 < 0) rel.r <- NaN
	
	#-------------------------------

	## average effects of both latent indicators
	eff2 <- (as.matrix(RR.1$effects[,-1]) + as.matrix(RR.2$effects[,-1])) / 2
	eff3 <- data.frame(RR.1$effects$id, eff2)
	colnames(eff3) <- colnames(RR.1$effects)
	eff <- eff3
	
	attr(eff[,2], "type") <- "actor"
	attr(eff[,3], "type") <- "partner"
	if (ncol(eff)>=4) attr(eff[,4], "type") <- "self"
	
	eff2 <- (as.matrix(RR.1$effects.gm[,-1]) + as.matrix(RR.2$effects.gm[,-1])) / 2
	eff3 <- data.frame(RR.1$effects.gm$id, eff2)
	colnames(eff3) <- colnames(RR.1$effects.gm)
	eff.gm <- eff3
	attr(eff.gm[,2], "type") <- "actor"
	attr(eff.gm[,3], "type") <- "partner"
	
	
	effRel2 <- merge(RR.1$effectsRel, RR.2$effectsRel, by=c("actor.id", "partner.id", "dyad"))
	effRel2$relationship <- apply(effRel2[,c("relationship.x", "relationship.y")], 1, mean, na.rm=TRUE)
	effRel <- effRel2[,c("actor.id", "partner.id", "dyad", "relationship")]
	effRel <- effRel[order(effRel$dyad),]
	
	
	attr(eff[,grep(localOptions$suffixes[1], colnames(eff), fixed=TRUE)], "reliability") <- rel.a
	attr(eff[,grep(localOptions$suffixes[2], colnames(eff), fixed=TRUE)], "reliability") <- rel.p
	attr(effRel$relationship, "reliability") <- rel.r
	
	res <- list(effects = eff, effects.gm=eff.gm, effectsRel=effRel, varComp=results, unstabdycov1=unstabdycov1, unstabper1=unstabper1, unstabtar1=unstabtar1, unstabrel1=unstabrel1, unstable.raw=unstable.raw, latent=TRUE, SEVAR=c(biSEVAR[1:2], biSEVAR["sesch2"], NA, (biSE["sesag"]^2+biSE["sesbf"]^2)/2, biSEVAR[6]), anal.type="Latent construct analysis of one construct measured by two round robin variables")
	attr(res, "group.size") <- n
	attr(res, "varname") <- paste(attr(RR.1, "varname"), attr(RR.2, "varname"), sep="/")
	if ((attr(RR.1, "self") == TRUE) & (attr(RR.2, "self") == TRUE)) {attr(res, "self") <- TRUE} else {attr(res, "self") <- FALSE}
	class(res) <- "RRuni"
}

	
	return(res)
}




# helper function
ifg <- function(g) {
	if (is.null(g)) {
		return("");
	} else {
		return(paste("|",g));
	}
}


# Wrapper function: depending on parameters, different results are calculated:
# @param se Either "SOREMO" (= between group significance test) or "LashleyBond" (= advanced significance test)
#' @export
#' @importFrom reshape2 melt
#' @importFrom reshape2 acast
#' @importFrom plyr ldply
#' @importFrom plyr laply
RR <- function(formula, data, na.rm=FALSE, minData = 1, verbose=TRUE, g.id=NULL, index="", exclude.ids="", varname=NA, se="LashleyBond", minVar=localOptions$minVar, ...) {

	extra <- list(...)
	
	se <- match.arg(se, choices=c("LashleyBond", "SOREMO"))

	# set default
	analysis <- "manifest"
	RRMatrix2 <- RRMatrix3 <- RRMatrix4 <- NULL

	# transform long format (formula) to quadratic matrices
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");

	#remove spaces from formula
	f1 <- formula
	lhs <- strsplit(gsub(" ","",as.character(f1)[2], fixed=TRUE), "+", fixed=TRUE)[[1]]
	rhs <- strsplit(gsub(" ","",as.character(f1)[3], fixed=TRUE),"\\*|\\|", perl=TRUE)[[1]]

	actor.id <- rhs[1]
	partner.id <- rhs[2]
	if (length(rhs)>=3) {group.id <- rhs[3]} else {group.id=NULL}


	# if a grouping factor is provided: forward to RR.multi
	if (!is.null(group.id)) {
		res <- RR.multi(f1, data=data, na.rm=na.rm, verbose=verbose, index=index, minData=minData, exclude.ids=exclude.ids, varname=varname, se=se, ...)
		
		if (!is.null(res$univariate)) {
			
			# bivariate case
			
			# if variance < minVar: set effects to NA
			if (!is.na(minVar)) {
				if (checkVar(res$univariate[[1]]$varComp[1, 3], minVar)) {
					res$univariate[[1]]$effects[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[1]]$effects.gm[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(res$univariate[[1]]$varComp[2, 3], minVar)) {
					res$univariate[[1]]$effects[,4][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[1]]$effects.gm[,4][1:nrow(res$univariate[[1]]$effects)] <- NA

				}
				if (checkVar(res$univariate[[2]]$varComp[1, 3], minVar)) {
					res$univariate[[2]]$effects[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[2]]$effects.gm[,3][1:nrow(res$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(res$univariate[[2]]$varComp[2, 3], minVar)) {
					res$univariate[[2]]$effects[,4][1:nrow(res$univariate[[1]]$effects)] <- NA
					res$univariate[[2]]$effects.gm[,4][1:nrow(res$univariate[[1]]$effects)] <- NA

				}
			}
			
			
		} else {
				# if variance < minVar: set effects to NA
				if (!is.na(minVar)) {
					if (checkVar(res$varComp[1, 3], minVar)) {
						res$effects[,3][1:nrow(res$effects)] <- NA
						res$effects.gm[,3][1:nrow(res$effects)] <- NA
					}
					if (checkVar(res$varComp[2, 3], minVar)) {
						res$effects[,4][1:nrow(res$effects)] <- NA
						res$effects.gm[,4][1:nrow(res$effects)] <- NA

					}
				}
		}
		
		res$minVar <- minVar
		for (g in 1:length(res$groups)) res$groups[[g]]$minVar <- minVar
		res$se <- se
		
		return(res)
		
	}


	# ---------------------------------------------------------------------
	#  Single group

	# univariater Fall:
	if (length(lhs)==1) {
		lhs1 <- strsplit(lhs, "/")[[1]]
	
		# manifester vs. latenter Fall
		if (length(lhs1)==1) {
			RRMatrix1 <- long2matrix(formula(paste(lhs1,"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=exclude.ids, ...)[[1]]
			analysis <- "manifest"
			
		} else if (length(lhs1)==2) {
			
			# if (!is.null(extra[["bistyle"]])) {v2 <- FALSE} else {v2=verbose}
			
			ex1 <- attr(long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex2 <- attr(long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex3 <- Reduce(union, list(exclude.ids, ex1, ex2))
			
			
			RRMatrix1 <- long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3, bistyle=TRUE)[[1]]
			RRMatrix2 <- long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3, ...)[[1]]
			analysis <- "latent"
		}
	
	} else 
	# bivariater Fall
	if (length(lhs)==2) {
	
		lhs1 <- strsplit(lhs[1], "/")[[1]]
	
		# manifester vs. latenter Fall		
		if (length(lhs1)==1) {
			# univariat:
			RRMatrix1 <- long2matrix(formula(paste(lhs[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=exclude.ids, ...)[[1]]
			
			RRMatrix2 <- long2matrix(formula(paste(lhs[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=exclude.ids, ...)[[1]]
			analysis <- "manifest"
		} else if (length(lhs1)==2) {
			# latent
			lhs2 <- strsplit(lhs[2], "/")[[1]]
			
			
			# exclude participants
			
			ex1a <- attr(long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex1b <- attr(long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex2a <- attr(long2matrix(formula(paste(lhs2[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex2b <- attr(long2matrix(formula(paste(lhs2[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, ...), "excluded.participants")
			ex3 <- Reduce(union, list(exclude.ids, ex1a, ex1b, ex2a, ex2b))
			
			RRMatrix1 <- long2matrix(formula(paste(lhs1[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=verbose, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3, bistyle=TRUE)[[1]]
			RRMatrix2 <- long2matrix(formula(paste(lhs1[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3)[[1]]
			RRMatrix3 <- long2matrix(formula(paste(lhs2[1],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3)[[1]]
			RRMatrix4 <- long2matrix(formula(paste(lhs2[2],"~",actor.id,"*",partner.id,ifg(group.id))), data, verbose=FALSE, minData=minData, skip3=TRUE, g.id=g.id, exclude.ids=ex3)[[1]]
			analysis <- "latent"
		}
	} else {stop("Error: Unknown term in formula.")}
	



## if all RRMatrices are NULL: stop
if (is.null(RRMatrix1) & is.null(RRMatrix2) & is.null(RRMatrix3) & is.null(RRMatrix4)) {
	return(NULL);
}
	
# depending on given parameters different results are calculated

#-----------------------------
#- One group

	if (is.null(RRMatrix2) & is.null(RRMatrix3) & is.null(RRMatrix4)) {
		if (analysis=="latent") {
			return(NULL);
			# warning("Warning: analysis='latent' only is valid, when two different RRMatrices for one latent construct are given")
		}
		
		res <- RR.univariate(RRMatrix1, na.rm, verbose, index=index, varname=varname)
		
		# if variance < minVar: set effects to NA
		if (!is.na(minVar)) {
			if (checkVar(res$varComp[1, 3], minVar)) {
				res$effects[,2][1:nrow(res$effects)] <- NA
				res$effects.gm[,2][1:nrow(res$effects)] <- NA
			}
			if (checkVar(res$varComp[2, 3], minVar)) {
				res$effects[,3][1:nrow(res$effects)] <- NA
				res$effects.gm[,3][1:nrow(res$effects)] <- NA
				
			}
		}
		
		res$minVar <- minVar
		res$se <- se
		return(res)
	}
	
#-----------------------------
#- Two groups, independent or latent constructs

	if (is.null(RRMatrix3) & is.null(RRMatrix4)) {
		
		if (is.null(RRMatrix1) | is.null(RRMatrix2)) {
			# if (verbose) {warning("Error: One of both round robin matrices has to few participants!", call.=FALSE)}
			return();
		}
		res <- RR.bivariate(RRMatrix1, RRMatrix2, analysis=analysis, na.rm=na.rm, verbose=verbose, index=index, varname=varname)
		
		
		if (!is.null(res$univariate)) {
			
			# bivariate case
			
			# if variance < minVar: set effects to NA
			if (!is.na(minVar)) {
				# Actor variance below minVar? Remove actor effects
				if (checkVar(res$univariate[[1]]$varComp[1, 3], minVar)) {
					res$univariate[[1]]$effects[,2] <- NA
					res$univariate[[1]]$effects.gm[,2] <- NA
				}
				# Partner variance below minVar? Remove partner effects
				if (checkVar(res$univariate[[1]]$varComp[2, 3], minVar)) {
					res$univariate[[1]]$effects[,3] <- NA
					res$univariate[[1]]$effects.gm[,3] <- NA

				}
				# Actor variance below minVar? Remove actor effects
				if (checkVar(res$univariate[[2]]$varComp[1, 3], minVar)) {
					res$univariate[[2]]$effects[,2] <- NA
					res$univariate[[2]]$effects.gm[,2] <- NA
				}
				# Partner variance below minVar? Remove partner effects
				if (checkVar(res$univariate[[2]]$varComp[2, 3], minVar)) {
					res$univariate[[2]]$effects[,3] <- NA
					res$univariate[[2]]$effects.gm[,3] <- NA

				}
			}
			
			
		} else {
				# if variance < minVar: set effects to NA
				if (!is.na(minVar)) {
					if (checkVar(res$varComp[1, 3], minVar)) {
						res$effects[,2] <- NA
						res$effects.gm[,2] <- NA
					}
					if (checkVar(res$varComp[2, 3], minVar)) {
						res$effects[,3] <- NA
						res$effects.gm[,3] <- NA
					}
				}
		}
		
		res$minVar <- minVar
		res$se <- se
		return(res);
	}
	
	
#-----------------------------
#- four groups: two constructs measured with each two variables

	if (!is.null(RRMatrix1) & !is.null(RRMatrix2) & !is.null(RRMatrix3) & !is.null(RRMatrix4)) {
		
		
		# calculate latent effects for both constructs
		lat1.full <- RR.bivariate(RRMatrix1, RRMatrix2, analysis="latent", na.rm=na.rm, verbose=FALSE, index=index, varname=varname)
		lat2.full <- RR.bivariate(RRMatrix3, RRMatrix4, analysis="latent", na.rm=na.rm, verbose=FALSE, index=index, varname=varname)
		lat1 <- lat1.full$varComp$estimate
		lat2 <- lat2.full$varComp$estimate
				
		
		# calculate new raw data: mean of both indicators; a manifest bivariate analysis is conducted on the mean variable
		# --> all calculations are correct, except the standardization --> this has to be done on the latent data
		RR12 <- (RRMatrix1+RRMatrix2)/2
		RR34 <- (RRMatrix3+RRMatrix4)/2
		
		#TAG1 <-- this is a bookmark, do not remove
		bivariate <- RR.bivariate(RR12, RR34, na.rm=na.rm, verbose=FALSE, noCorrection=TRUE)$bivariate
		
		#Estimation of bivariate relations on construct level
		
		w <- getOption("warn")
		options(warn=-1)
		denom <- c(
			sqrt(lat1[1]*lat2[1]),
			sqrt(lat1[2]*lat2[2]),
			sqrt(lat1[1]*lat2[2]),
			sqrt(lat1[2]*lat2[1]),
			sqrt (lat1[3]*lat2[3]),
			sqrt (lat1[3]*lat2[3])
		)
		options(warn=w)
		
		bivariate$standardized <- clamp(bivariate$estimate / denom)

		# erase covariances if one variance component is < 0
		if (lat1[1] <= 0) bivariate[c(1,3), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (lat1[2] <= 0) bivariate[c(2,4), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (lat2[1] <= 0) bivariate[c(1,4), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (lat2[2] <= 0) bivariate[c(2,3), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (lat1[3] <= 0) bivariate[c(5,6), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (lat2[3] <= 0) bivariate[c(5,6), c("standardized", "se", "t.value", "p.value")] <- NaN
		
		
		univariate <- list()
		univariate[[1]] <- lat1.full
		univariate[[2]] <- lat2.full
		
		grandres <- list(univariate = univariate, bivariate = bivariate)
		class(grandres) <- "RR"
		grandres$anal.type <- "Bivariate analysis of two constructs, each measured by two round robin variables"
		attr(grandres, "group.size") <- nrow(RRMatrix2)
		
		
			# if variance < minVar: set effects to NA
			if (!is.na(minVar)) {
				if (checkVar(grandres$univariate[[1]]$varComp[1, 3], minVar)) {
					grandres$univariate[[1]]$effects[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[1]]$effects.gm[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(grandres$univariate[[1]]$varComp[2, 3], minVar)) {
					grandres$univariate[[1]]$effects[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[1]]$effects.gm[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA

				}
				if (checkVar(grandres$univariate[[2]]$varComp[1, 3], minVar)) {
					grandres$univariate[[2]]$effects[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[2]]$effects.gm[,3][1:nrow(grandres$univariate[[1]]$effects)] <- NA
				}
				if (checkVar(grandres$univariate[[2]]$varComp[2, 3], minVar)) {
					grandres$univariate[[2]]$effects[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA
					grandres$univariate[[2]]$effects.gm[,4][1:nrow(grandres$univariate[[1]]$effects)] <- NA

				}
			}

		grandres$minVar <- minVar
		grandres$se <- se
		return(grandres)		
	} else {
		# warning("Error: One of the round robin matrices has to few participants!", call.=FALSE)
	}
	
	return(NULL);
}


 


# uni1, uni2: univariate Analysen der beiden Konstrukte (Daten werden zum Standardisieren der bivariaten Koeffizienten gebraucht)
getWTest <- function(RR0, res1, typ="univariate", uni1=NA, uni2=NA, unstable=NA, se="LashleyBond") {
	
	if (is.null(RR0)) return();
	
	if (typ=="univariate") {
		if (length(RR0$univariate)==2) {varComp <- RR0$univariate[[1]]$varComp} else {varComp <- RR0$varComp}
		
		# create empty template for varComp
		varComp[, -1] <- NA

		for (v in names(table(res1$type))) {
			
			#------------------------------------
			#--  Significance test for Variance Components (OLD STYLE)
			# -- calculate weighted mean and weighted between groups t-test
			#------------------------------------
			
			if (se == "SOREMO") {
				w.t <- weighted.t.test(res1$estimate[res1$type == v], res1$group.size[res1$type == v]-1, mu=0)
				varComp[varComp$type==v,]$estimate <- w.t$estimate
				varComp[varComp$type==v,]$se <- w.t$se
				varComp[varComp$type==v,]$SEVAR <- w.t$se^2
				varComp[varComp$type==v,]$t.value <- w.t$statistic
				varComp[varComp$type==v,]$p.value <- w.t$p.value
			}
			
			#------------------------------------
			#--  Significance test for Variance Components (LASHLEY-BOND STYLE)
			#------------------------------------
			
			if (se == "LashleyBond") {
				SEs2 <- res1$SEVAR[res1$type == v]
				VAR <- res1$estimate[res1$type == v]
				
				# compute weights based on n
				w <- res1$group.size[res1$type == v]-1

				# weighted mean estimate of the SRM parameter				
				VAR.mean.weighted <- sum(VAR*w) / sum(w)
				
				SE.mean.weighted <- sqrt(sum(w^2*SEs2) / (sum(w)^2))
					
				t.value <- VAR.mean.weighted / SE.mean.weighted
								
				df <- sum(res1$group.size[res1$type == v]-1)
				# compute two-tailed p-value (p-values for variances are divided by two later)
 				p.value <- pt(abs(t.value), df, lower.tail=FALSE)*2
				
				varComp[varComp$type==v,]$estimate <- VAR.mean.weighted
				varComp[varComp$type==v,]$se <- as.vector(SE.mean.weighted)
				varComp[varComp$type==v,]$SEVAR <- as.vector(SE.mean.weighted)^2
				varComp[varComp$type==v,]$t.value <- t.value
				varComp[varComp$type==v,]$p.value <- p.value
			}
		}


		varComp$p.value[1:4] <- varComp$p.value[1:4] / 2
		# Varianzen nur einseitig testen (Voreinstellung bei weighted.t.test ist zweiseitig)
		
		# unstable variance im latent bivariaten Fall wird von aussen in die Funktion gegeben
		if (!is.na(unstable)) {
			varComp[4,2] <- unstable
		}

		rownames(varComp) <- NULL

		#standardized coefficients need special treatment ...
		varComp[1,3] <- posOrNA(varComp[1,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		varComp[2,3] <- posOrNA(varComp[2,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		varComp[3,3] <- posOrNA(varComp[3,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		varComp[4,3] <- posOrNA(varComp[4,2])/ sum(posOrNA(varComp[1:4,2]), na.rm=TRUE)
		w <- getOption("warn")
		options(warn=-1)
			varComp[5,3] <- varComp[5,2]/ sqrt(varComp[1,2]*varComp[2,2])
			varComp[6,3] <- varComp[6,2]/ varComp[3,2]
		options(warn=w)
		
		varComp[,3] <- clamp(varComp[,3])
		
		# variance below zero: erase all other indices
		bz <- which(varComp[1:3, 2]<0)
		if (length(bz)>0) {
			varComp[bz, c("standardized", "se", "t.value", "p.value")] <- NaN
			if (varComp[1,2]<0 | varComp[2,2]<0) varComp[5, c("standardized", "se", "t.value", "p.value")] <- NaN
		}
		
		# error variance: set se, t, p to NA (instead of NaN)
		varComp[4,4:6] <- NA
	
		return(varComp)
	}
	
	if (typ=="bivariate") {
		bivariate <- RR0$bivariate
		bivariate[, -1] <- NA
		#bivariate$p.value <- NA

		for (v in names(table(res1$type))) {
			#------------------------------------
			#--  Significance test for Variance Components (OLD STYLE)
			# -- calculate weighted mean and weighted between groups t-test
			#------------------------------------
			
			if (se == "SOREMO") {
				w.t <- weighted.t.test(res1$estimate[res1$type == v], res1$group.size[res1$type == v]-1, mu=0)
				bivariate[bivariate$type==v,]$estimate <- w.t$estimate
				bivariate[bivariate$type==v,]$se <- w.t$se
				bivariate[bivariate$type==v,]$biSEVAR <- w.t$se^2
				bivariate[bivariate$type==v,]$t.value <- w.t$statistic
				bivariate[bivariate$type==v,]$p.value <- w.t$p.value
			}
			
			#------------------------------------
			#--  Significance test for Variance Components (LASHLEY-BOND STYLE)
			#------------------------------------
			
			if (se == "LashleyBond") {
				SEs2 <- res1$biSEVAR[res1$type == v]
				VAR <- res1$estimate[res1$type == v]
				
				# compute weights based on n
				w <- res1$group.size[res1$type == v]-1

				# weighted mean estimate of the SRM parameter				
				VAR.mean.weighted <- sum(VAR*w) / sum(w)
				
				SE.mean.weighted <- sqrt(sum(w^2*SEs2) / (sum(w)^2))
					
				t.value <- VAR.mean.weighted / SE.mean.weighted
								
				df <- sum(res1$group.size[res1$type == v]-1)
 				p.value <- pt(t.value, df, lower.tail=FALSE)
				
				bivariate[bivariate$type==v,]$estimate <- VAR.mean.weighted
				bivariate[bivariate$type==v,]$se <- as.vector(SE.mean.weighted)
				bivariate[bivariate$type==v,]$biSEVAR <- as.vector(SE.mean.weighted)^2
				bivariate[bivariate$type==v,]$t.value <- t.value
				bivariate[bivariate$type==v,]$p.value <- p.value
			}			
		}
		
		#standardized coefficients need special treatment ...
		w <- getOption("warn")
		options(warn=-1)
			bivariate[1,3] <- bivariate[1,2]/ sqrt(uni1[1,2]*uni2[1,2])
			bivariate[2,3] <- bivariate[2,2]/ sqrt(uni1[2,2]*uni2[2,2])
			bivariate[3,3] <- bivariate[3,2]/ sqrt(uni1[1,2]*uni2[2,2])
			bivariate[4,3] <- bivariate[4,2]/ sqrt(uni1[2,2]*uni2[1,2])
			bivariate[5,3] <- bivariate[5,2]/ sqrt(uni1[3,2]*uni2[3,2])
			bivariate[6,3] <- bivariate[6,2]/ sqrt(uni1[3,2]*uni2[3,2])
		options(warn=w)
		
		
		# erase covariances if one variance component is < 0
		if (uni1[1,2] <= 0) bivariate[c(1,3), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (uni2[1,2] <= 0) bivariate[c(1,4), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (uni1[2,2] <= 0) bivariate[c(2,4), c("standardized", "se", "t.value", "p.value")] <- NaN
		if (uni2[2,2] <= 0) bivariate[c(2,3), c("standardized", "se", "t.value", "p.value")] <- NaN

		bivariate[,3] <- clamp(bivariate[,3])
	
		return(bivariate)
	}
}




RR.multi.uni <- function(formule, data, na.rm=FALSE, verbose=TRUE, index="", minData=1, exclude.ids="", varname=NA, se="LashleyBond", ...) {

	# this function needs data in long format ...
	extra <- list(...)
		
	# parse formula
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");
	
	# f1 = formula without grouping factor
	fstring <- paste(as.character(formule)[c(2,1,3)], collapse=" ")
	f0 <- strsplit(gsub(" ","",fstring, fixed=TRUE),"\\|", perl=TRUE)[[1]]
	f1 <- formula(f0[1])
	f3 <- strsplit(strsplit(gsub(" ","",fstring, fixed=TRUE),"~", perl=TRUE)[[1]][1], "+", fixed=TRUE)[[1]]
	group.id <- f0[2]

	mode <- ifelse(length(f3)==2,"bi","uni")
		
	res <- data.frame()
	res.bi <- data.frame()
	g.uni <- list()
	
	# n.m stores the group sizes
	saa <- sbb <- scc <- sccs <- n.m <- c()
	sesaa2 <- sesbb2 <- sescc2 <- sesab2 <- sesccs2 <- c()	
	undc1 <- unp1 <- unt1 <- unr1 <- un.raw  <- c()
	
	self <- FALSE	# are self ratings present?
	
	for (g in names(table(data[,group.id]))) {
		
		#print(g)
		RR0 <- RR(f1, data=data[data[,group.id] == g,], verbose=verbose, na.rm=na.rm, g.id=group.id, index=index, minData=minData, exclude.ids=exclude.ids, varname=varname, minVar=NA, ...)
		
		#print(str(RR0))

		if (is.null(RR0)) {next;} else {RR1 <- RR0}
		if (attr(RR0, "self") == TRUE) {self <- TRUE}
		g.id <- g
		
		RR0$group.id <- g.id

		g.uni[[g]] <- RR0
		
		if (RR0$latent==FALSE) {
			saa <- c(saa, RR0$varComp[1,2])
			sbb <- c(sbb, RR0$varComp[2,2])
			scc <- c(scc, RR0$varComp[3,2])
			sccs <- c(sccs, RR0$varComp[6,2])
		} else {
			undc1 <- c(undc1, RR0$unstabdycov1)
			unp1 <- c(unp1, RR0$unstabper1)
			unt1 <- c(unt1, RR0$unstabtar1)
			unr1 <- c(unr1, RR0$unstabrel1)
		}
		
		n.m <- c(n.m, attr(RR0, "group.size"))
		
		u1 <- RR0$varComp
		u1$SEVAR <- RR0$SEVAR
		u1$variable <- 1

		u1$group.size <-  attr(RR0, "group.size")
		u1$group.id <- g.id
		
		res <- rbind(res, u1)
		
	}
	
	#stop()
	
	# aus der liste die Effekte extrahieren und zusammenfuegen
	effect <- ldply(g.uni, function(x) {return(x$effects)})
	effect[,1:2] <- effect[,2:1]
	colnames(effect)[1:2] <- c("id", "group.id")
	effect[,1] <- factor(effect[,1])
	effect[,2] <- factor(effect[,2])
	
	type <- c("actor", "partner", "self")
	for (ty in 3:ncol(effect)) {
		attr(effect[,ty], "type") <- type[ty-2]
	}

	eff.gm <- ldply(g.uni, function(x) {return(x$effects.gm)})
	eff.gm[,1:2] <- eff.gm[,2:1]
	colnames(eff.gm)[1:2] <- c("id", "group.id")
	eff.gm[,1] <- factor(eff.gm[,1])
	eff.gm[,2] <- factor(eff.gm[,2])
	
	effectRel <- ldply(g.uni, function(x) {return(x$effectsRel)})
	colnames(effectRel)[1:3] <- c("group.id", all.vars(f1)[2:3])
	
	effectRel[,1] <- factor(effectRel[,1])
	effectRel[,2] <- factor(effectRel[,2])
	effectRel[,3] <- factor(effectRel[,3])
	effectRel[,4] <- factor(effectRel[,4])

	# im latenten Fall: die Error variance erst am Ende berechnen (d.h., alle error componenten ueber alle Gruppen mitteln, die unter NUll auf Null setzen, dann addieren)
	
	unstable.raw.m <- NA
	if (RR1$latent==TRUE) {
		unstable.raw.m <- max(weighted.mean(unp1, n.m), 0) + max(weighted.mean(unt1, n.m), 0) + max(weighted.mean(unr1, n.m), 0)
	}	

	if (length(effect) == 0) {
		effect <- data.frame(actor=NA, partner=NA, relationship=NA)
	}
	
	# get weighted variance components
	varComp <- getWTest(RR1, res, unstable=ifelse(is.null(unstable.raw.m), NULL, unstable.raw.m), se=se)

	# calculate reliability for actor and partner effects, and variance of group means
	group.var <- NA
	
	if (!is.null(n.m)) {

		n <- mean(n.m)
		
		if (RR1$latent==FALSE) {
			saa.m <- weighted.mean(saa, n.m-1)
			sbb.m <- weighted.mean(sbb, n.m-1)
			scc.m <- weighted.mean(scc, n.m-1)
			sccs.m <- weighted.mean(sccs, n.m-1)
	
			rel.a <- saa.m / (saa.m + scc.m*(n-1)/(n*(n-2)) + sccs.m/(n*(n-2)))
			if (saa.m < 0) rel.a <- NaN
			rel.p <- sbb.m / (sbb.m + scc.m*(n-1)/(n*(n-2)) + sccs.m/(n*(n-2)))
			if (sbb.m < 0) rel.p <- NaN
		} else {
			
			
			unp1.m <- weighted.mean(unp1, n.m-1, na.rm=TRUE)
			unt1.m <- weighted.mean(unt1, n.m-1, na.rm=TRUE)
			unr1.m <- weighted.mean(unr1, n.m-1, na.rm=TRUE)
			undc1.m <- weighted.mean(undc1, n.m-1, na.rm=TRUE)
			
			
			r <- 2 # r = number of replications - in our case, it's always 2
			rel.a <- varComp$estimate[1] / ((varComp$estimate[1] + (unp1.m/r)) + (varComp$estimate[3]+(unr1.m/r))*(n-1)/(n*(n-2)) + (varComp$estimate[6]+(undc1.m/r))/(n*(n-2)))
			if (varComp$estimate[1] < 0) rel.a <- NaN

			
			rel.p <- varComp$estimate[2] / ((varComp$estimate[2] + (unt1.m/r)) + (varComp$estimate[3]+(unr1.m/r))*(n-1)/(n*(n-2)) + (varComp$estimate[6]+(undc1.m/r))/(n*(n-2)))
			if (varComp$estimate[2] < 0) rel.p <- NaN

			rel.r <- varComp$estimate[3] / (varComp$estimate[3]+(unr1.m/r))
			if (varComp$estimate[3] < 0) rel.r <- NaN
			
			attr(effectRel$relationship, "reliability") <- clamp(rel.r)
		}

		attr(effect[,grep(localOptions$suffixes[1], colnames(effect), fixed=TRUE)], "reliability") <- clamp(rel.a)
		attr(effect[,grep(localOptions$suffixes[2], colnames(effect), fixed=TRUE)], "reliability") <- clamp(rel.p)
		
		# group mean & group variance		
		group.means <- tapply(data[, all.vars(f1)[1]], data[, group.id], mean, na.rm=TRUE)
		group.var <- var(group.means) - ((varComp$estimate[1] + varComp$estimate[2] + 2*varComp$estimate[5])/n + (varComp$estimate[3]+varComp$estimate[6])/(n*(n-1)))
				
	}	
	
	
	if (se == "LashleyBond") {
		ST <- "(significance test based on Lashley & Bond, 1997, Psychological Methods)"
	}
	if (se == "SOREMO") {
		ST <- "(significance test based on SOREMO; Kenny & LaVoie, 1984)"
	}
	anal.type <- paste0(RR1$anal.type, " in multiple groups ", ST)

	if (!is.null(varComp)) {

		#res2 <- list(effects = effect, effectsRel = effectRel, effects.gm = eff.gm, varComp = varComp, groups = g.uni, varComp.groups=res, group.var=group.var, biSEVAR=varComp$se^2, anal.type=anal.type)
		res2 <- list(effects = effect, effectsRel = effectRel, effects.gm = eff.gm, varComp = varComp, groups = g.uni, varComp.groups=res, group.var=group.var, anal.type=anal.type)
		class(res2) <- "RRmulti"
		attr(res2, "varname") <- attr(g.uni[[1]], "varname")
		attr(res2, "self") <- self
	
		# # noch rausfinden, welche Teilnehmer ausgeschlossen wurden
		# l1 <- long2matrix(formule, data, verbose=FALSE)
		# attr(res2, "excluded.participants") <- attr(l1, "excluded.participants")
		# attr(res2, "excluded.groups") <- attr(l1, "excluded.groups")
	
		return(res2)
	} else {return();}
}




# @param se Either "SOREMO" (= between group significance test) or "LashleyBond" (= advanced significance test)
RR.multi <- function(formule, data, na.rm=FALSE, verbose=TRUE, index="", minData=1, exclude.ids="", varname=NA, se="LashleyBond", ...) {

	# this function needs data in long format ...
	extra <- list(...)
	
	# parse formula
	if (is.null(data)) stop("If a formula is specified, an explicit data object has to be provided!");
	
	# f1 = formula without grouping factor
	fstring <- paste(as.character(formule)[c(2,1,3)], collapse=" ")
	f0 <- strsplit(gsub(" ","",fstring, fixed=TRUE),"\\|", perl=TRUE)[[1]]
	f1 <- formula(f0[1])
	group.id <- f0[2]
	
	f3 <- strsplit(strsplit(gsub(" ","",fstring, fixed=TRUE),"~", perl=TRUE)[[1]][1], "+", fixed=TRUE)[[1]]
	f4 <- strsplit(gsub(" ","",fstring, fixed=TRUE),"~", perl=TRUE)[[1]][2]
	
	if (sum(grepl("/", f3, fixed=TRUE))>1) {analysis <- "latent"} else {analysis <- "manifest"}

	# Vorlage fuer die Varianzkomponente erstellen
	df <- data[data[,group.id]==data[1,group.id],]
	mode <- ifelse(length(f3)==2,"bi","uni")

	if (mode=="uni") {
		return(RR.multi.uni(formule, data, na.rm, verbose, index=index, minData=minData, exclude.ids=exclude.ids, varname=varname, se=se, ...))
	}



	# ... ansonsten bi-mode durchfuehren
	
	# find out, which participants are excluded and exclude them from all variables
	if (analysis=="manifest") {
		ex1 <- attr(long2matrix(formula(paste(f3[1], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex2 <- attr(long2matrix(formula(paste(f3[2], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex3 <- Reduce(union, list(ex1, ex2, exclude.ids))
	} else {
		# latent analysis
		ex1a <- attr(long2matrix(formula(paste(strsplit(f3[1], "/", fixed=TRUE)[[1]][1], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex1b <- attr(long2matrix(formula(paste(strsplit(f3[1], "/", fixed=TRUE)[[1]][2], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex2a <- attr(long2matrix(formula(paste(strsplit(f3[2], "/", fixed=TRUE)[[1]][1], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex2b <- attr(long2matrix(formula(paste(strsplit(f3[2], "/", fixed=TRUE)[[1]][2], "~", f4)), data, verbose=FALSE, minData=minData), "excluded.participants")
		ex3 <- Reduce(union, list(ex1a, ex1b, ex2a, ex2b, exclude.ids))
	}
	
	V1 <- RR.multi.uni(formula(paste(f3[1], "~", f4)), data, na.rm, verbose=verbose, index=index, minData=minData, exclude.ids=ex3, bistyle=TRUE, se=se)
	V2 <- RR.multi.uni(formula(paste(f3[2], "~", f4)), data, na.rm, verbose=FALSE, index=index, minData=minData, exclude.ids=ex3, se=se)
	V2$varComp.groups$variable <- 2
		
	res.bi <- data.frame()
	bi.groups <- list()
	
	for (g in names(table(data[,group.id]))) {
		
			RR0 <- RR(f1, data=data[data[,group.id] == g,], verbose=FALSE, na.rm=na.rm, minData=minData, exclude.ids=ex3, minVar=NA, se=se)
			
			if (is.null(RR0)) {next;} else
			{RR1 <- bi.groups[[g]]  <- RR0}
			
			if (!is.null(RR1$bivariate)) {
				res.bi <- rbind(res.bi, data.frame(RR1$bivariate, group.size=attr(RR1, "group.size"), group=g))
			}
		
	}
	
	if (se == "LashleyBond") {
		ST <- "(significance test based on Lashley & Bond, 1997, Psychological Methods)"
	}
	if (se == "SOREMO") {
		ST <- "(significance test based on SOREMO; Kenny & LaVoie, 1984)"
	}
	anal.type <- paste0(RR1$anal.type, " in multiple groups ", ST)
	
	bivariate <- getWTest(RR1, res.bi, typ="bivariate", V1$varComp, V2$varComp, se=se)
	
	res <- list(univariate = list(V1, V2), bivariate = bivariate, SEVAR=bivariate$se^2, anal.type=anal.type, groups=bi.groups)
	class(res) <- "RRmulti"

	return(res)
}

