############################
## Compute distance to center for a group
############################

dissassoc <- function(diss, group, weights=NULL, R=1000, weight.permutation="replicate", squared=FALSE){
	return(dissassocweighted(diss, group, weights, R, weight.permutation, squared))
}
olddissassoc <- function(diss, group , R=1000) {
	#Notation comme pour l'ANOVA, SC=Inertia dans le sens du criète de Ward
	if (inherits(diss, "dist")) {
		diss <- dist2matrix(diss)
	}
	dissmatrix <- diss[!is.na(group), !is.na(group)]
	SCtot <- .Call(TMR_tmrsubmatrixinertia, dissmatrix, 1:nrow(dissmatrix))
	n <- nrow(dissmatrix)
	ind <- 1:n
	grp <- factor(group[!is.na(group)])
	lgrp <- levels(grp)
	k <- length(lgrp)
	#pour chaque valeur du groupe
	#ret$contrib <- vector("list", length(lgrp))
	indgrp <- list()
	ret <- list()
	ret$groups <- data.frame(n=numeric(k+1), variance=numeric(k+1))
	rownames(ret$groups) <- c(lgrp, "Total")
	SCres <- 0
	for (i in 1:k) {
		#on crée le groupe en question
		cond <- grp==lgrp[i]
		ret$groups$n[i] <- sum(cond)
		if (ret$groups$n[i]==0) {
			indgrp[i] <- integer()
			ret$groups$variance[i] <- 0
		}
		else {
			indgrp[[i]] <- ind[cond]
			r <- .Call(TMR_tmrsubmatrixinertia, dissmatrix, as.integer(indgrp[[i]]))
			ret$groups$variance[i] <- r/ret$groups$n[i]
			SCres <- r+SCres
		}
	}

	ret$groups$variance[length(lgrp)+1] <- SCtot/nrow(dissmatrix)
	ret$groups$n[length(lgrp)+1] <- nrow(dissmatrix)
	ret$anova.table <- data.frame(SS=c(SCtot-SCres, SCres, SCtot),
			df=c(k-1, n-k, n-1), MSE=c((SCtot-SCres)/(k-1),
			SCres/(n-k), SCtot/(n-1)))
	rownames(ret$anova.table) <- c("Exp", "Res", "Total")
#	print(indgrp)
#	print(SCtot)
	if (R==1) {
		vals <- internalBootstrapCompareGroups(ind, ind, dissmatrix, indgrp, SCtot)
		ret$stat <- data.frame(PseudoF=vals[1], PseudoR2=vals[2],
				PseudoF_Pval=NA, PseudoT=vals[3], PseudoT_Pval=NA)
		ret$perms <- NA
	}
	else {
		bts <- boot(ind, internalBootstrapCompareGroups, R, sim="permutation", stype="i",
			dissmatrix=dissmatrix, indgrp=indgrp, SCtot=SCtot)
		ret$stat <- data.frame(PseudoF=bts$t0[1], PseudoR2=bts$t0[2],
			PseudoF_Pval=sum(bts$t[, 1]>bts$t0[1])/bts$R,
			PseudoT=bts$t0[3], PseudoT_Pval=sum(bts$t[, 3]>bts$t0[3])/bts$R)
		ret$perms <- bts
	}
	rownames(ret$stat) <- c("")
	ret$call <- match.call()
	ret$R <- R
	class(ret) <- "olddissassoc"
	return(ret)
}
print.olddissassoc <- function(x, ...) {
	cat("Pseudo ANOVA table:\n")
	print(x$anova.table, ...)
	cat("\nTest values ", "(p-values based on", (x$R-1), "permutations):\n")
	print(x$stat, ...)
	cat("\nVariance per level:\n")
	print(x$groups, ...)
}

hist.olddissassoc <- function(x, test="PseudoF", breaks="FD", main=paste("Distribution of", test), xlab=test, pvalue.limit=NULL, freq=FALSE, ...) {
	if (test=="PseudoF") {
		ti <- 1
	}
	else if (test=="PseudoT") {
		ti <- 3
	}
	else {
		stop("test argument should be one of PseudoF or PseudoT")
	}
	if (x$R==1) {
		stop("Cannot plot permutation test distribution for R = 1")
	}
	testbootorder <- order(x$perms$t[, ti], decreasing=TRUE)
	hist(x$perms$t[, ti], main=main, xlab=xlab, breaks=breaks, freq=freq, ...)
	if (!is.null(pvalue.limit)) {
		abline(v=x$perms$t[testbootorder[round(pvalue.limit*x$R)], ti], col="blue")
	}
	abline(v=x$stat[, test], col="red")
}


internalBootstrapCompareGroups <- function(seqdata, ind, dissmatrix, indgrp, SCtot) {
	SCres <- 0
	#pour chaque valeur du groupe
	k <- length(indgrp)
	n <- length(seqdata)
	s1ni <- 0
	lns <- 0
	nlnvi <- 0
	for (i in 1:k) {
		#on crée le groupe en question
		groupe <- sort.int(ind[as.integer(indgrp[[i]])], method="quick")
#		SCresi <- .Call("tmrsubmatrixinertia", dissmatrix, groupe, PACKAGE="TraMineR")
#		groupe <- as.integer(ind[indgrp[[i]]])
		ni <- length(groupe)
		#on calcul l'inertie intraclasse
		SCresi <- .Call(TMR_tmrsubmatrixinertia, dissmatrix, groupe)
		vari <- SCresi/ni
		s1ni <- s1ni+(1/(ni-1))
		lns <- lns+(ni-1)*(vari/(n-k))
		nlnvi <- nlnvi+(ni-1)*log(vari)
		SCres <- SCres+SCresi
	}
	Tcalc <- (n-k)*log(lns)-nlnvi
	Ccalc <- 1+1/(3*(k-1))*(s1ni-1/(n-k))
	Bartlett <- Tcalc/Ccalc
#		retourne le gain d'inertie
	SCexp <- SCtot-SCres
	PseudoR2 <- SCexp/SCtot
	PseudoF <- (SCexp/(k-1))/(SCres/(n-k))
	return(c(PseudoF, PseudoR2, Bartlett))
}
