############################
## Compute distance to center for a group
############################
disscenter <- function(diss, group=NULL, medoids.index=NULL, allcenter=FALSE, weights=NULL, squared=FALSE) {
disscentertrim(diss=diss, group=group, medoids.index=medoids.index, allcenter=allcenter, weights=weights, squared=squared, trim =0)
}
disscentertrim <- function(diss, group=NULL, medoids.index=NULL, allcenter=FALSE, weights=NULL, squared=FALSE, trim =0) {
	if(is.logical(medoids.index)){
		if(medoids.index){
			medoids.index <- "First"
		}else{
			medoids.index <- NULL
		}
	}
	retmedoids <- !is.null(medoids.index)
	if (retmedoids) {
		allcenter <- FALSE
	}
	allmedoids <- FALSE
	if(!is.null(medoids.index)) {
		if(medoids.index=="all"){
			allmedoids <- TRUE
		} else if (medoids.index!="first") {
			stop('medoids.index argument should be one of "first", "all" or NULL')
		}
	}
	
	
	## max.iter <- 20
	if (inherits(diss, "dist")) {
		diss <- as.matrix(diss)
	}
	if (squared) {
		diss <- diss^2
	}
	if (is.null(group)) {
		group <- integer(nrow(diss))
		group[] <- 1
	}
	ind <- 1:nrow(diss)
	grp <- factor(group)
	lgrp <- levels(grp)
	if (is.null(weights)) {
		weights <- rep(1, nrow(diss))
	}
	weights <- as.double(weights)
	if(allcenter){
		ret <- data.frame(numeric(nrow(diss)))
	}else{
		ret <- numeric(nrow(diss))
	}
	if (retmedoids) {
		if (allmedoids) {
			medoids <- list()
		}
		else {
			medoids <- numeric(length(lgrp))
		}
	}
	keep=1-trim
	## pour chaque valeur du groupe
	for (i in 1:length(lgrp)) {
		## on crée le groupe en question
		cond <- grp==lgrp[i]
		grpindiv <- sort(ind[cond])
		## on calcul la contribution a l'inertie intraclasse
		if (allcenter) {
			ret[, i] <- 0
			others <- sort(ind[!cond])
			dT <- .Call(TMR_tmrWeightedInertiaContribExt, diss, grpindiv, others,weights)
			dTindiv <- 1:sum(cond)
			dT <- dT - weighted.mean(dT[dTindiv], weights[grpindiv])
			ret[grpindiv, i] <- dT[dTindiv]
			ret[others, i] <- dT[-(dTindiv)]
		}
		else {
			dc <- .Call(TMR_tmrWeightedInertiaContrib, diss, as.integer(grpindiv),weights)
			dc <- dc-weighted.mean(dc, weights[cond])/2
			
			ret[grpindiv] <- dc
			mindc <- min(dc)
			if (trim>0) {
					maxdist <- quantile(dc, probs=keep)
					trimmedcond <- dc<=maxdist
					dT <- .Call(TMR_tmrinertiacontribext, diss, grpindiv[trimmedcond], grpindiv[!trimmedcond])
					ntrimmed <- sum(trimmedcond)
					dT <- dT-mean(dT[1:ntrimmed])/2
					mindc <- min(dT)
					ret[grpindiv[trimmedcond]] <- dT[1:ntrimmed]
					ret[grpindiv[!trimmedcond]] <- dT[(ntrimmed+1):length(grpindiv)]
			}
			if (retmedoids) {
				if (allmedoids) {
					medoids[[i]] <- which(ret==mindc & cond)
				}
				else {
					medoids[i] <- sort(which(ret==mindc & cond))[1]
				}
			}
		}
	}
	if (retmedoids) {
		## No group
		if (length(lgrp)==1) {
			return(medoids[[1]])
		}
		names(medoids) <- lgrp
		return(medoids)
	}
	if (allcenter) {
		names(ret) <- lgrp
	}
	return(ret)
}