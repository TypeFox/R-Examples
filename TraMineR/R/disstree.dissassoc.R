############################
## Compute distance to center for a group
############################


DTNdissassocweighted <- function(dmat, grp, indiv, R, weights, weight.permutation) {
	dmatsize <-  as.integer(nrow(dmat))
	internalDTNdissassocUnweighted <- function(grp, ind, indiv, dmat, grp2, use.sort) {
		groupe1 <- indiv[ind[grp]]
		groupe2 <- indiv[ind[grp2]]
		 if (use.sort) {
				groupe1 <- sort.int(groupe1, method="quick")
				groupe2 <- sort.int(groupe2, method="quick")
		 }
		 r1 <- .Call(TMR_tmrsubmatrixinertia, dmat, as.integer(groupe1))
		 r2 <- .Call(TMR_tmrsubmatrixinertia, dmat, as.integer(groupe2))
		 return(-(r1+r2))
	}
	internalDTNdissassocWeighted <- function(grp, ind, indiv, dmat, grp2, use.sort, weights, permutGroup) {
		if (permutGroup) {
			weights[indiv] <- weights[indiv[ind]]
		}
		groupe1 <- indiv[ind[grp]]
		groupe2 <- indiv[ind[grp2]]
		 if (use.sort) {
				groupe1 <- sort.int(groupe1, method="quick")
				groupe2 <- sort.int(groupe2, method="quick")
		 }
		 r1 <- .Call(TMR_tmrWeightedInertiaDist, dmat, dmatsize,
				as.integer(FALSE), as.integer(groupe1), as.double(weights), 
				as.integer(FALSE))
		r2 <- .Call(TMR_tmrWeightedInertiaDist, dmat, dmatsize,
				as.integer(FALSE), as.integer(groupe2), as.double(weights), 
				as.integer(FALSE))
		 return(-(r1+r2))
	}
	internalDTNdissassocReplicate <- function(grp, ind, indiv, dmat, grp2) {
		groupe1n <- indiv[ind[grp]]
		groupe2n <- indiv[ind[grp2]]
		wwt1 <- tabulate(groupe1n)
		wwt2 <- tabulate(groupe2n)
		groupe1<- which(wwt1>0)
		groupe2<- which(wwt2>0)
		r1 <- .Call(TMR_tmrWeightedInertiaDist, dmat, dmatsize,
				as.integer(FALSE), as.integer(groupe1), as.double(wwt1), 
				as.integer(FALSE))
		r2 <- .Call(TMR_tmrWeightedInertiaDist, dmat, dmatsize, 
				as.integer(FALSE), as.integer(groupe2), as.double(wwt2), 
				as.integer(FALSE))
		return(-(r1+r2))
	}
	if ( weight.permutation %in% c("diss", "group")) {
		perms <- TraMineR.permutation(grp, R=R, statistic=internalDTNdissassocWeighted, 
			dmat=dmat, indiv=indiv, grp2=(!grp), use.sort=(length(grp)>750), weights=weights, 
			permutGroup=(weight.permutation=="group"))
	} else if(weight.permutation=="none") {
		perms <- TraMineR.permutation(grp, R=R, statistic=internalDTNdissassocUnweighted, 
			dmat=dmat, indiv=indiv, grp2=(!grp), use.sort=(length(grp)>750))
	}else { ##Replicate
		grp <- rep(grp, times=as.integer(weights[indiv]))
		indiv <- rep(indiv, times=as.integer(weights[indiv]))
		perms <- TraMineR.permutation(grp, R=R, statistic=internalDTNdissassocReplicate, 
			dmat=dmat, indiv=indiv, grp2=(!grp))
	}
	return(perms$pval[1])
}


