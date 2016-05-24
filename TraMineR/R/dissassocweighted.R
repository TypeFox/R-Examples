############################
## Compute distance to center for a group
############################
WeightedResidualsZ <- function(valg, weightedvalg){
	return(sum(weightedvalg*(valg-weighted.mean(valg, weightedvalg))^2))
}
UnweightedResidualsZ <- function(valg){
	return(sum((valg-mean(valg))^2))
}
ComputeWeightedTestValues <- function(allindiv, ind, SCresiFunc, weights, SCtot, totweights, k, disscSCtot){
		SCres <- 0
		#pour chaque valeur du groupe
		lns <- 0
		nlnvi <- 0
		FBFdenomin <- 0
		sumz <- 0
		s1ni <- 0
		#sumexp <- 0
		for (i in 1:k) {
			#on crée le groupe en question
			SCresif <- SCresiFunc(allindiv, ind, i, weights)
			SCresi <- SCresif[2]
			ni <- SCresif[1]
			##sumexp <- sumexp+ SCresif[3]
			sumz <- sumz+SCresif[3]
			## SCresi <- .Call("tmrsubmatrixinertia", dissmatrix, groupe, PACKAGE="TraMineR")
			vari <- SCresi/ni
			lns <- lns+(ni-1)*(vari/(totweights-k))
			nlnvi <- nlnvi+(ni-1)*log(vari)
			SCres <- SCres+SCresi
			FBFdenomin <- FBFdenomin+ (1-ni/totweights)*vari
			s1ni <- s1ni+(1/(ni-1))
		}
		Ccalc <- 1+1/(3*(k-1))*(s1ni-1/(totweights-k))
		Tcalc <- (totweights-k)*log(lns)-nlnvi
		Bartlett <- Tcalc/Ccalc
	#		retourne le gain d'inertie
		SCexp <- SCtot-SCres
		PseudoR2 <- SCexp/SCtot
		PseudoF <- (SCexp/(k-1))/(SCres/(totweights-k))
		FPF<-SCexp/FBFdenomin
		PseudoW <- ((disscSCtot-sumz)/(k-1))/(sumz/(totweights-k))
		#print(c(PseudoF, PseudoR2, Bartlett))
		return(c(PseudoF, FPF, PseudoR2, Bartlett, PseudoW))
	}

TraMineR.indgrplist <- function(grpint, k=max(grpint)){
	## building group information once, before permutation to save computation
	## indgrp contain indice of indrep of each group
	indgrp <- list()
	## Compute for each group the variances table
	for (i in 1:k) {
		# Group belong condition
		cond <- grpint==i
		if (sum(cond)>0) {
			indgrp[[i]] <- which(cond)
		}
	}
	return(indgrp)
}

dissassocweighted.permgroup <- function(diss, grpint, weights, R, ret, randomWeight) {
	# SCtot <- ret$anova.table$SS[3]
	k <- length(ret$groups$n)-1
	totweights <- ret$groups$n[k+1]
	indgrp <- TraMineR.indgrplist(grpint, k)
	n <- nrow(diss)
	localstuffpermgroup <- new.env()
	localstuffpermgroup$computeT0 <- TRUE
	WeightedGroupTestValues <- function(allindiv, ind){
		##permut weights
		if(randomWeight) {
			if(localstuffpermgroup$computeT0) {
				localstuffpermgroup$computeT0 <- FALSE
				ww <- weights
			}
			else {
				ww <- weights[sample(1:n)]
			}
		} else{
			ww <- weights[ind]
		}
		
		SCtot <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), 
			as.integer(FALSE), allindiv, as.double(ww), 
			as.integer(FALSE))
		dissc <- disscenter(diss, group=grpint[ind], weights=ww)
		disscSCtot <- WeightedResidualsZ(dissc, ww)
		
		SCresiNoReplicate <- function(allindiv, ind, i, ww){
			wind <- sort.int(ind[as.integer(indgrp[[i]])], method="quick")
			ni <- sum(ww[as.integer(indgrp[[i]])])
			SCresi <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), 
				as.integer(FALSE), as.integer(wind), as.double(ww), 
				as.integer(FALSE))
			wind <- indgrp[[i]]
			wwi <- ww[wind]
			return(c(ni, SCresi, WeightedResidualsZ(dissc[wind], wwi)))
		}
		return(ComputeWeightedTestValues(allindiv=allindiv, ind=ind, 
			SCresiFunc=SCresiNoReplicate, weights=ww, SCtot=SCtot, totweights=totweights, 
			k=k, disscSCtot=disscSCtot))
		
	}
	return(TraMineR.permutation(1:n, R, WeightedGroupTestValues))
	
}
dissassocweighted.permdiss <- function(diss, grpint, weights, R, ret) {
	SCtot <- ret$anova.table$SS[3]
	k <- length(ret$groups$n)-1
	totweights <- ret$groups$n[k+1]
	indgrp <- TraMineR.indgrplist(grpint, k)
	n <- nrow(diss)
	dissc <- disscenter(diss, group=grpint, weights=weights)
	disscSCtot <- WeightedResidualsZ(dissc, weights)

	WeightedGroupTestValues <- function(allindiv, ind){
		SCresiNoReplicate <- function(allindiv, ind, i, ww){
			wind <- sort.int(ind[as.integer(indgrp[[i]])], method="quick")
			ni <- sum(ww[as.integer(indgrp[[i]])])
			SCresi <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), 
				as.integer(FALSE), as.integer(wind), as.double(ww), 
				as.integer(FALSE))
			wwi <- ww[wind]
			return(c(ni, SCresi, WeightedResidualsZ(dissc[wind], wwi)))
		}
		return(ComputeWeightedTestValues(allindiv=allindiv, ind=ind, 
			SCresiFunc=SCresiNoReplicate, weights=weights, SCtot=SCtot, totweights=totweights, 
			k=k, disscSCtot=disscSCtot))
		
	}
	return(TraMineR.permutation(1:n, R, WeightedGroupTestValues))
	
}

dissassocweighted.unweighted <- function(diss, grpint, weights, R, ret) {
	SCtot <- ret$anova.table$SS[3]
	k <- length(ret$groups$n)-1
	totweights <- ret$groups$n[k+1]
	indgrp <- TraMineR.indgrplist(grpint, k)
	n <- nrow(diss)
	dissc <- disscenter(diss, group=grpint)
	disscSCtot <- UnweightedResidualsZ(dissc)

	WeightedGroupTestValues <- function(allindiv, ind){
		SCresiUnweighted <- function(allindiv, ind, i, ww){
			##all indiv is indrep
			wind <- sort.int(ind[as.integer(indgrp[[i]])], method="quick")
			ni <- length(wind)
			SCresi <- .Call(TMR_tmrsubmatrixinertia, diss, 
				wind)
			return(c(ni, SCresi, UnweightedResidualsZ( dissc[wind])))
		}
		return(ComputeWeightedTestValues(allindiv=allindiv, ind=ind, 
			SCresiFunc=SCresiUnweighted, weights=weights, SCtot=SCtot, totweights=totweights, 
			k=k, disscSCtot=disscSCtot))
		
	}
	return(TraMineR.permutation(1:n, R, WeightedGroupTestValues))
	
}

dissassocweighted.replicate <- function(diss, grpint, weights, R, ret) {
	SCtot <- ret$anova.table$SS[3]
	k <- length(ret$groups$n)-1
	totweights <- ret$groups$n[k+1]
	n <- nrow(diss)
	indrep <- as.integer(rep(1:n, times=as.integer(weights)))
	grouprep=as.integer(rep(grpint, times=as.integer(weights)))
	indgrp <- TraMineR.indgrplist(grouprep, k)
	dissc <- disscenter(diss, group=grpint, weights=weights)
	disscSCtot <- WeightedResidualsZ(dissc, weights)
	

	WeightedGroupTestValues <- function(allindiv, ind){
		SCresiReplicate <- function(allindiv, ind, i, ww){
			##all indiv is indrep
			nind <- allindiv[ind[as.integer(indgrp[[i]])]]
			wwt <- tabulate(nind)
			wind <- which(wwt>0)
			ni <- length(nind)
			SCresi <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), 
				as.integer(FALSE), as.integer(wind), as.double(wwt), 
				as.integer(FALSE))
			wwti <- wwt[wind]
			return(c(ni, SCresi, WeightedResidualsZ(dissc[wind], wwti)))
		}
		return(ComputeWeightedTestValues(allindiv=allindiv, ind=ind, 
			SCresiFunc=SCresiReplicate, weights=weights, SCtot=SCtot, totweights=totweights, 
			k=k, disscSCtot=disscSCtot))
		
	}
	return(TraMineR.permutation(indrep, R, WeightedGroupTestValues))
}

dissassocweighted.internaldisscenter <- function(diss, grpint, indiv, k, weights){
	alldc <- numeric(length(indiv))
	for(i in 1:k) {
		cond <- grpint==i
		if(sum(cond)>0){
			grpindiv <- indiv[cond]
			dc <- .Call(TMR_tmrWeightedInertiaContrib, diss, as.integer(grpindiv), weights)
			alldc[grpindiv] <- dc-weighted.mean(dc, weights[grpindiv])/2
		}
	}
	return(alldc)
}

dissassocweighted.permbootstrap <- function(diss, grpint, weights, R, ret, samplesize, sampleprob) {
	# SCtot <- ret$anova.table$SS[3]
	k <- length(ret$groups$n)-1
	totweights <- ret$groups$n[k+1]
	n <- nrow(diss)
	wwall <- rep(1, n)
	# allpop <- 1:samplesize
	WeightedGroupTestValues <- function(dissindiv, grpindiv){
		##permut weights
		grpp <- grpint[grpindiv]
		SCtot <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), 
			as.integer(FALSE), dissindiv, as.double(wwall), 
			as.integer(FALSE))
		dissc <- dissassocweighted.internaldisscenter(diss, grpp, dissindiv,k, wwall)
		disscSCtot <- UnweightedResidualsZ(dissc)
		
		SCresiNoReplicate <- function(allindiv, ind, i, ww){
			indgrp <- grpp==i
			wind <- sort.int(dissindiv[indgrp], method="quick")
			ni <- length(wind)
			SCresi <- .Call(TMR_tmrWeightedInertiaDist, diss, as.integer(n), 
				as.integer(FALSE), as.integer(wind), as.double(wwall), 
				as.integer(FALSE))
			#SCresi <- .Call("tmrWeightedInertiaDist", dissmatrix, as.integer(n), 
			#	as.integer(FALSE), as.integer(wind), as.double(ww), 
			#	as.integer(FALSE), PACKAGE="TraMineR")
			return(c(ni, SCresi, UnweightedResidualsZ(dissc[indgrp])))
		}
		return(ComputeWeightedTestValues(allindiv=dissindiv, ind=grpindiv, 
			SCresiFunc=SCresiNoReplicate, weights=wwall, SCtot=SCtot, totweights=totweights, 
			k=k, disscSCtot=disscSCtot))
		
	}
	return(TraMineR.permutationweight(grpint, R=R, statistic=WeightedGroupTestValues,samplesize=samplesize, sampleprob=sampleprob,
	t0=dissassocweighted.permdiss(diss=diss, grpint=grpint, weights=weights, R=0, ret=ret)$t0))
	
}

dissassocweighted <- function(diss, group, weights, R, weight.permutation, squared, samplesize=NULL) {
	
	## Notation comme pour l'ANOVA, SC=Inertia dans le sens du criète de Ward
	if (inherits(diss, "dist")) {
		diss <- dist2matrix(diss)
	}
	## Manage missing values by removing them
	dissmatrix <- diss[!is.na(group), !is.na(group)]
	
	n <- nrow(dissmatrix)
	if(is.null(samplesize)) {
		samplesize <- n
	}
	if (squared) {
		dissmatrix <- dissmatrix^2
	}
	
	use.replicate <- weight.permutation %in% c("rounded-replicate", "replicate")
	unweighted <- is.null(weights)

	if (unweighted) {
		weights <- rep(1, n)
		use.replicate <- FALSE
		weight.permutation <- "none"
	} else {
		weights <- as.double(weights[!is.na(group)])	
	}
	## Allow integer weights for replicates
	if(use.replicate) {
		rounderror <- sum(abs(round(weights, 0) - weights))
		if(rounderror>0){
			if (weight.permutation=="replicate") {
				stop(" [!] to permute replicate, you should specify integer weights")
			}
			message("Weigths loss : ", rounderror, " (", (rounderror/sum(weights)), ")")
			weights <- round(weights, 0)
		}
		weights <- as.integer(weights)
	}
	if(weight.permutation=="random-sampling"){
		w2 <- (weights*samplesize)/sum(weights)
		rounderror <- mean(abs(w2 - weights))
		if(rounderror>0){
			message(" [>] weights corrected to match sample size. Mean correction ", format(rounderror))
		}
		weights <- w2
	}
	grp <- factor(group[!is.na(group)])
	
	grpint <- as.integer(grp)
	## Number of groups
	k <- length(levels(grp))
	allindiv <- 1:n
	## Compute basic global statistics used after
	totweights <- sum(weights)
	SCtot <- .Call(TMR_tmrWeightedInertiaDist, dissmatrix, as.integer(n), 
		as.integer(FALSE), allindiv, as.double(weights), 
		as.integer(FALSE))
	
	## Building return the returned object
	ret <- list()
	ret$groups <- data.frame(n=numeric(k+1), discrepancy=numeric(k+1))
	rownames(ret$groups) <- c(levels(grp), "Total")
	
	## Residual sum of square
	SCres <- 0
	## Compute for each group the discrepancy table
	for (i in 1:k) {
		# Group belong condition
		cond <- grpint==i
		## Group size
		ret$groups$n[i] <- sum(weights[cond])
		## Empty group ?
		if (ret$groups$n[i]==0) {
			ret$groups$discrepancy[i] <- 0
		}else {
			## Intra group residual
			r <- .Call(TMR_tmrWeightedInertiaDist, dissmatrix, as.integer(n), 
				as.integer(FALSE), as.integer(which(cond)), as.double(weights), 
				as.integer(FALSE))
			ret$groups$discrepancy[i] <- r/ret$groups$n[i]
			SCres <- r+SCres
		}
	}
	## Total discrepancy
	ret$groups$discrepancy[k+1] <- SCtot/totweights
	ret$groups$n[k+1] <- totweights
	ret$anova.table <- data.frame(SS=c(SCtot-SCres, SCres, SCtot),
			df=c(k-1, totweights-k, totweights-1), MSE=c((SCtot-SCres)/(k-1),
			SCres/(totweights-k), SCtot/(totweights-1)))
	rownames(ret$anova.table) <- c("Exp", "Res", "Total")
	#print(indgrp)
	## Starting permutation test
	## Computing approximated SCtot by expansion
	## For Bartlett test
	if(unweighted) {
		ret$perms <- dissassocweighted.unweighted(diss=dissmatrix, grpint, weights=weights, R=R, ret=ret)
	}
	else if (use.replicate){
		ret$perms <- dissassocweighted.replicate(diss=dissmatrix, grpint, weights=weights, R=R, ret=ret)
	}
	else if(weight.permutation=="diss") {
		ret$perms <- dissassocweighted.permdiss(diss=dissmatrix, grpint, weights=weights, R=R, ret=ret)
	}
	else if(weight.permutation=="group"){
		ret$perms <- dissassocweighted.permgroup(diss=dissmatrix, grpint, weights=weights, R=R, ret=ret, randomWeight=FALSE)
	}
	else if(weight.permutation=="permute"){
		ret$perms <- dissassocweighted.permgroup(diss=dissmatrix, grpint, weights=weights, R=R, ret=ret, randomWeight=TRUE)
	}
	else if(weight.permutation=="random-sampling"){
		ret$perms <- dissassocweighted.permbootstrap(diss=dissmatrix, grpint, weights=weights, R=R, ret=ret, samplesize=samplesize, sampleprob=weights/sum(weights))
	}
	
	
	ret$stat <- summary(ret$perms)	
	ret$weight.permutation <- weight.permutation
	rownames(ret$stat) <- c("Pseudo F", "Pseudo Fbf", "Pseudo R2", "Bartlett", "Levene")
	ret$call <- match.call()
	ret$R <- R
	class(ret) <- "dissassoc"
	return(ret)
}
print.dissassoc <- function(x, pvalue.confint=0.95, digits = NULL, ...) {
	if(is.null(digits)) {
		digits = getOption("digits")
	}
	cat("Pseudo ANOVA table:\n")
	print(x$anova.table, digits=digits, ...)
	cat("\nTest values ", "(p-values based on", (x$R), "permutation):\n")
	if (x$weight.permutation!="none"){
		cat("Weights permuted using \"", x$weight.permutation, "\"\n")
	}
	stat <- x$stat
	print(stat,digits=digits, ...)
	if(!is.null(pvalue.confint)){
		## Two sided
		cat("\nInconclusive intervals: \n")
		pvalue.confint <- (1-(1-pvalue.confint)/2)
		pvalues <- c(0.01, 0.05)
		confinter <- qnorm(pvalue.confint)*sqrt(pvalues*(1-pvalues)/x$R)
		cat(paste(format(pmax(0,(pvalues-confinter)), digits=3), " < ",pvalues, " < ",format(pmin(1,(pvalues+confinter)),digits=3), collapse="\n"), "\n")
	}
	
	cat("\nDiscrepancy per level:\n")
	print(x$groups, digits=digits, ...)
}

hist.dissassoc <- function(x, test="Pseudo F", ...) {
	ti <- match(test, rownames(x$stat))
	hist(x$perms, index=ti, xlab=test,...)
}


