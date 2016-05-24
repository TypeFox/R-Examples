dltMatchCurvePoints <- function(lm.list, cal.coeff, min.direct.tangency = 25, 
	min.fill.tangency = 0, epi.err.weight = 1, rec.err.weight = 0){

	# CHECK IF LIST INPUT IS MORE THAN ONE LANDMARK/CURVE
	if(!is.null(names(lm.list))){

		# EMPTY LISTS FOR STORING OUTPUT VALUES
		epipolar_dist <- list()
		curve_pt_dist <- list()
		portion_matched_nonref <- list()

		for(landmark_name in names(lm.list)){

			# CALL MATCH CURVE POINTS FOR EACH LANDMARK/CURVE SET
			dlt_match_curve_points <- dltMatchCurvePoints(lm.list=lm.list[[landmark_name]], 
				cal.coeff=cal.coeff, min.direct.tangency=min.direct.tangency,  
				min.fill.tangency=min.fill.tangency, epi.err.weight=epi.err.weight, 
				rec.err.weight=rec.err.weight)

			# SAVE OUTPUT VALUES TO LISTS
			lm.list[[landmark_name]] <- dlt_match_curve_points$match.lm.list
			epipolar_dist[[landmark_name]] <- dlt_match_curve_points$epipolar.dist
			match_count[[landmark_name]] <- dlt_match_curve_points$match.count
		}
		
		r <- list(match.lm.list=lm.list, epipolar.dist=epipolar_dist, match.count=match_count)
		class(r) <- 'dltMatchCurvePoints'
		return(r)
	}

	# CHECK THAT THERE ARE AT LEAST TWO CAMERA VIEWS
	if(length(lm.list) < 2){r <- list(match.lm.list=lm.list);class(r) <- 'dltMatchCurvePoints';return(r)}		

	# CHECK THAT THERE ARE MORE THAN TWO POINTS IN EVERY VIEW
	pt_count <- rep(NA, length(lm.list))
	for(i in 1:length(lm.list)){
		if(is.matrix(lm.list[[i]])){
			pt_count[i] <- nrow(lm.list[[i]])
			if(nrow(lm.list[[i]]) < 3){r <- list(match.lm.list=lm.list);class(r) <- 'dltMatchCurvePoints';return(r)}
		}else{
			r <- list(match.lm.list=lm.list);class(r) <- 'dltMatchCurvePoints';return(r)
		}
	}

	# FIND MAX OR MIN FOR REFERENCE VIEW IF SPECIFIED
	ref.view <- which.max(pt_count)
	
	# SWITCH COLUMNS IF FIRST VIEW IS NOT REFERENCE VIEW
	if(ref.view == 2){
		lm.list <- list(lm.list[[2]], lm.list[[1]])
		cal.coeff <- cal.coeff[, ncol(cal.coeff):1]
	}

	# HOMLOGOUS CURVE POINT LISTS
	homol_curve_pts_list <- list()
	homol_curve_pts_list[[1]] <- matrix(lm.list[[1]], nrow=nrow(lm.list[[1]]), ncol=ncol(lm.list[[1]]), dimnames=list(rownames(lm.list[[1]]), colnames(lm.list[[1]])))
	homol_curve_pts_list[[2]] <- matrix(NA, nrow=nrow(lm.list[[1]]), ncol=ncol(lm.list[[1]]), dimnames=list(rownames(lm.list[[1]]), colnames(lm.list[[1]])))

	# ASSUME START AND END POINTS ARE HOMOLOGOUS
	homol_curve_pts_list[[2]][1, ] <- lm.list[[2]][1, ]
	homol_curve_pts_list[[2]][nrow(lm.list[[1]]), ] <- lm.list[[2]][nrow(lm.list[[2]]), ]

	# MATCH CURVE POINTS
	res <- matchCurvePoints(lm.list[[1]], lm.list[[2]], cal.coeff, min.direct.tangency = min.direct.tangency, 
		min.fill.tangency = min.fill.tangency, epi.err.weight = epi.err.weight, rec.err.weight = rec.err.weight)

	# SAVE RESULTING MATRICES
	#homol_curve_pts_list[[1]] <- lm.list[[1]][!is.na(res$match[, 1]), ]
	homol_curve_pts_list[[2]] <- res$match

	# GET EPIPOLAR DISTANCES BETWEEN REFERENCE AND TEST POINTS
	epipolar_dist <- dltEpipolarDistance(p1=homol_curve_pts_list[[1]], p2=homol_curve_pts_list[[2]], cal.coeff, reciprocal=TRUE)
	
	# SWITCH BACK POINT LISTS IF FIRST VIEW IS NOT REFERENCE VIEW
	if(ref.view == 2) homol_curve_pts_list <- list(homol_curve_pts_list[[2]], homol_curve_pts_list[[1]])

	r <- list(
		match.lm.list=homol_curve_pts_list, 
		epipolar.dist=epipolar_dist, 
		match.count=res$match.cts
	)
	class(r) <- 'dltMatchCurvePoints'

	return(r)
}

summary.dltMatchCurvePoints <- function(object, print.tab = '', ...){

	r <- '\n'
	r <- c(r, print.tab, 'dltMatchCurvePoints Summary\n')

	if(!is.null(names(object$match.lm.list))){
		curve_found <- F

		for(curve_name in names(object$match.lm.list)){
			if(is.null(object$epipolar.dist[[curve_name]])) next
			r <- c(r, print.tab, '\tCurve name: ', curve_name, '\n')
			#r <- c(r, print.tab, '\tReference point count: ', length(object$epipolar.dist[[curve_name]]), '\n')
			r <- c(r, print.tab, '\tNumber matched initial/fill/unmatched: ', paste(object$match.count[[curve_name]], collapse='/'), '\n')
			r <- c(r, print.tab, '\tStart/end epipolar distances: ', round(object$epipolar.dist[[curve_name]][1], 2), ' px, ', round(object$epipolar.dist[[curve_name]][length(object$epipolar.dist[[curve_name]])], 2), ' px\n')
			r <- c(r, print.tab, '\tEpipolar distances, Mean: ', round(mean(object$epipolar.dist[[curve_name]], na.rm=TRUE), 2), ' px; Min: ', round(min(object$epipolar.dist[[curve_name]], na.rm=TRUE), 2), ' px; Max: ', round(max(object$epipolar.dist[[curve_name]], na.rm=TRUE), 2), ' px\n')
#			r <- c(r, print.tab, '\t\tPortion of non-reference matched: ', round(object$portion.matched.nonref[[curve_name]]*100, 1), ' %\n')
#			r <- c(r, print.tab, '\t\tMean epipolar-curve point distance: ', round(mean(object$curve.pt.dist[[curve_name]], na.rm = TRUE), 2), ' px +/- ', round(sd(object$curve.pt.dist[[curve_name]], na.rm = TRUE), 2), '\n')
#			r <- c(r, print.tab, '\t\tMax epipolar-curve point distance: ', round(max(object$curve.pt.dist[[curve_name]], na.rm = TRUE), 2), ' px\n')
			curve_found <- T
		}
		if(!curve_found) r <- c(r, print.tab, '\tNo curves found.\n')
	}else{
		if(length(object$epipolar.dist) == 0){
			r <- c(r, print.tab, '\tNo curves found.\n')
			class(r) <- "summary.dltMatchCurvePoints"
			return(r)
		}
	
		#r <- c(r, print.tab, '\tReference point count: ', length(object$epipolar.dist), '\n')
		r <- c(r, print.tab, '\tNumber matched initial/fill/unmatched: ', paste(object$match.count, collapse='/'), '\n')
		r <- c(r, print.tab, '\tStart/end epipolar distances: ', round(object$epipolar.dist[1], 2), ' px, ', round(object$epipolar.dist[length(object$epipolar.dist)], 2), ' px\n')
		r <- c(r, print.tab, '\tEpipolar distances, Mean: ', round(mean(object$epipolar.dist, na.rm=TRUE), 2), ' px; Min: ', round(min(object$epipolar.dist, na.rm=TRUE), 2), ' px; Max: ', round(max(object$epipolar.dist, na.rm=TRUE), 2), ' px\n')
#		r <- c(r, print.tab, '\t\tPortion of non-reference matched: ', round(object$portion.matched.nonref*100, 1), ' %\n')
#		r <- c(r, print.tab, '\t\tMean epipolar-curve point distance: ', round(mean(object$curve.pt.dist, na.rm = TRUE), 2), ' px +/- ', round(sd(object$curve.pt.dist, na.rm = TRUE), 2), '\n')
#		r <- c(r, print.tab, '\t\tMax epipolar-curve point distance: ', round(max(object$curve.pt.dist, na.rm = TRUE), 2), ' px\n')
	}

	class(r) <- "summary.dltMatchCurvePoints"
	r
}

print.summary.dltMatchCurvePoints <- function(x, ...) cat(x, sep='')
