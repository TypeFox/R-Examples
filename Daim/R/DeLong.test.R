
############################################################
# x: Matrix(nrow=length(labels)); 
# ref: Index d. Referenz-Markers (bei paarweisen Vergleichen; 
# für globalen p-Wert ist Referenz egal)


deLong.test <- function(x, labels, labpos, ref = NULL, conf.level = 0.95) 
{
	if(length(labels) != dim(x)[1])
		stop("\n The number of rows in x must match the length of labels\n")
	id.pos <- labels == labpos
	if(sum(id.pos) < 1)
		stop("\n wrong level specified!\n")
	if(dim(x)[2] < 2)
		stop("\n x must contain at least two columns!\n")
	if(dim(x)[1] < 2)
		stop("\n x must contain at least two rows!\n")
	nn <- sum(!id.pos)
	np <- sum(id.pos)
	nauc <- ncol(x)
	
	if(is.null(ref)) { 
		L <- matrix(0, nrow=nauc*(nauc-1)/2, ncol=nauc)
		newa <- 0
		for(i in 1:(nauc-1)) {
			newl <- nauc - i
			L[(newa+1):(newa+newl),i] <- rep(1, newl)
			L[(newa+1):(newa+newl),((i+1):(i+newl))] <- diag(-1, nrow=newl, ncol=newl)
			newa <- newa + newl
		}
	}
	else { 
		# test for superiority of one method against all others)
		if(ref > nauc) 
			stop(paste("Reference ref must be one of the markers (1...", nauc, " in this case)", sep=""))
		L <- matrix(1, ncol=nauc, nrow=nauc-1)
		L[,-ref] <- diag(-1, nrow=nauc-1, ncol=nauc-1)
	}

	markern <- as.matrix(x[!id.pos,])
	markerp <- as.matrix(x[id.pos,])

	###
	### compute wilcox statistic
	###
	WK.STAT <- function(x,y){
        r <- rank(c(x, y))
        n.x <- length(x)
        n.y <- length(y)
        STATISTIC <- sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2
        STATISTIC
	}
	
	auc <- vector("numeric", length=nauc)
	for(r in 1:nauc) {
		auc[r] <- WK.STAT(markerp[,r], markern[,r])
	}
	auc <- auc/(nn*np)

	###
	### if AUCs smaller than 0.5: 1-auc
	###
	if(any(auc < 0.5)) {
		x[,auc < 0.5] <- -x[,auc<0.5]
		auc[auc < 0.5] <- 1 - auc[auc<0.5]
		markern <- as.matrix(x[!id.pos,])
		markerp <- as.matrix(x[id.pos,])
	}

	V10 <- matrix(0, nrow=np, ncol=nauc)
	V01 <- matrix(0, nrow=nn, ncol=nauc)
	
	tmn <- t(markern)
	tmp <- t(markerp)
	for(i in 1:np) {
		V10[i,] <- rowSums(tmn < tmp[,i]) + 0.5 * rowSums(tmn == tmp[,i])
	}
	for(i in 1:nn) {
		V01[i,] <- rowSums(tmp > tmn[,i]) + 0.5 * rowSums(tmp == tmn[,i])
	}
	V10 <- V10/nn
	V01 <- V01/np

	W10 <- cov(V10)
	W01 <- cov(V01)
	
	###
	### estimated covariance matrix
	###
	S <- W10/np + W01/nn
	
	###
	### compute variances of AUCs and test for AUC > 0.5
	###
	
	### Hanley, McNeil (1982)
	q1 <- auc / (2 - auc)
	q2 <- 2*auc^2 / (1 + auc)
	
	### Haney, McNeil (1982) / Bamber (1975)
	aucvar <- (auc*(1 - auc) + (np - 1)*(q1 - auc^2) + (nn - 1)*(q2 - auc^2)) / (np*nn)
	zhalf <- (auc - 0.5) / sqrt(aucvar)
	phalf <- 1 - pnorm(zhalf)
	zdelong <- (auc - 0.5) / sqrt(diag(S))
	pdelong <- 1 - pnorm(zdelong)
	
	
	### global p-value
	###
	aucdiff <- L %*% auc
	z <- t(aucdiff) %*% matinv(L %*% S %*% t(L)) %*% aucdiff
	p <- pchisq(z, df=qr(L %*% S %*% t(L))$rank, lower.tail=FALSE)

	if(is.null(ref)) {
		cor.auc <- matrix(ncol=1, nrow=nauc*(nauc-1)/2)
		ci <- matrix(ncol=2, nrow=nauc*(nauc-1)/2)
		ctr <- 1
		rows <- vector("character", length=(nauc*(nauc-1)/2))
		pairp <- matrix(nrow=nauc*(nauc-1)/2, ncol=1)
		quantil <- qnorm(1 - (1 - conf.level)/2)
		for(i in 1:(nauc-1)) {
			for(j in (i+1):nauc) {
				cor.auc[ctr] <- S[i,j] / sqrt(S[i,i]*S[j,j])
				LSL <- t(c(1,-1)) %*% S[c(j,i),c(j,i)] %*% c(1,-1)
				tmpz <- (aucdiff[ctr]) %*% matinv(LSL) %*% aucdiff[ctr]
				pairp[ctr] <- 1 - pchisq(tmpz, df=qr(LSL)$rank)
				ci[ctr,] <- c(aucdiff[ctr] - quantil*sqrt(LSL), aucdiff[ctr] + quantil*sqrt(LSL))
				rows[ctr] <- paste(i, j, sep=" vs. ")
				ctr <- ctr+1
			}
		}
	} else {
		cor.auc <- matrix(ncol=1, nrow=nauc-1)
		ci <- matrix(ncol=2, nrow=nauc-1)
		rows <- vector("character", length=nauc-1)
		pairp <- matrix(nrow=nauc-1, ncol=1)
		comp <- (1:nauc)[-ref]
		for(i in 1:(nauc-1)) {
			cor.auc[i] <- S[ref,comp[i]] / sqrt(S[ref,ref] * S[comp[i],comp[i]])
			LSL <- t(c(1,-1)) %*% S[c(ref,comp[i]),c(ref,comp[i])] %*% c(1,-1)
			tmpz <- aucdiff[i] %*% matinv(LSL) %*% aucdiff[i]
			pairp[i] <- 1 - pchisq(tmpz, df=qr(LSL)$rank)
			ci[i,] <- c(aucdiff[i] - quantil*sqrt(LSL), aucdiff[i] + quantil*sqrt(LSL))
			rows[i] <- paste(ref, comp[i], sep=" vs. ")
		}		
	}

	newres <- as.data.frame(cbind(aucdiff, ci, pairp, cor.auc))
	names(newres) <- c("AUC Difference", "CI(lower)", "CI(upper)", "P.Value", "Correlation")
	rownames(newres) <- rows
	row.names(ci) <- row.names(cor.auc) <- row.names(aucdiff) <- row.names(pairp) <- rows
	colnames(ci) <- c(paste0(100*conf.level, "% CI (lower)"), paste0(100*conf.level, "% CI (upper)"))
	names(auc) <- 1:nauc
	auc <- as.data.frame(cbind(auc, sqrt(aucvar), phalf, sqrt(diag(S)), pdelong))
	colnames(auc) <- c("AUC", "SD(Hanley)", "P(H0: AUC=0.5)", "SD(DeLong)", "P(H0: AUC=0.5)")

	ERG <- list(AUC = auc, difference = newres, covariance = S, global.z = z, global.p = p)
	class(ERG) <- "DeLong"
	ERG
}


