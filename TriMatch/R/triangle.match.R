#' Creates matched triplets.
#' 
#' Create matched triplets by minimizing the total distance between matched triplets
#' within a specified caliper.
#' 
#' The \code{\link{trips}} function will estimate the propensity scores
#' for three models. This method will then find the best matched triplets based
#' upon minimizing the summed differences between propensity scores across the
#' three models. That is, the algorithm works as follows:
#' 
#' \itemize{
#' \item The first subject from model 1 is selected.
#' \item The \code{nmatch[1]} smallest distances are selected using propensity scores from
#'  model 1.
#' \item For each of the matches identified, the subjects propensity score from model
#'  2 is retrieved.
#' \item The \code{nmatch[2]} smallest distances are selected using propensity score from
#'  model 3.
#' \item For each of those matches identified, the subjects propensity score from model
#'  2 is retrieved.
#' \item The distances is calculated from the first and last subjects propensity scores
#'  from model 2.
#' \item The three distances are summed.
#' \item The triplet with the smallest overall distance is selected and returned.
#' }
#' 
#' @param tpsa the results from \code{\link{trips}}
#' @param caliper a vector of length one or three indicating the caliper to use 
#'        for matching within each step. This is expressed in standardized units such 
#'        that .25 means that matches must be within .25 of one standard deviation 
#'        to be kept, otherwise the match is dropped.
#' @param nmatch number of closest matches to retain before moving to next edge. This can
#'        be \code{Inf} in which case all matches within the caliper will be retained
#'        through to the next step. For large datasets, evaluating all possible
#'        matches within the caliper could be time consuming. 
#' @param match.order character vector of length three indicating the order in 
#'        which the matching algorithm will processes. The default is to use start
#'        with the group the middle number of subjects, followed by the smallest,
#'        and then the largest. 
#' @param exact a vector or data frame of representing covariates for exact matching.
#'        That is, matched triplets will first be matched exactly on these covariates
#'        before evalutating distances.
#' @param method This is a function that specifies which matched triplets will be
#'        retained. If \code{NULL}, all matched triplets within the specified
#'        caliper will be returned (equivelent to caliper matching in two group
#'        matching). The default is \code{\link{maximumTreat}} that
#'        attempts include each treatment at least once. 
#'        Another option is \code{\link{OneToN}} which mimicks the one-to-n 
#'        matching where treatments are matched to multiple control units.
#' @param ... other parameters passed to \code{method}.
#' @export
#' @examples
#' \dontrun{
#' data(turoing)
#' formu <- ~ Gender + Ethnicity + Military + ESL + EdMother + EdFather + Age +
#'      Employment + Income + Transfer + GPA
#' tpsa <- trips(tutoring, tutoring$treat, formu)
#' tmatch <- trimatch(tpsa, status=FALSE)
#' }
trimatch <- function(tpsa, 
					 caliper=.25, 
					 nmatch=c(15),
					 match.order,
					 exact,
					 method=maximumTreat, 
					 ...) {
	if(length(nmatch) == 1) {
		nmatch <- c(nmatch, nmatch)
	}
	if(length(caliper) == 1) {
		caliper <- rep(caliper, 3)
	}
	
	groups <- attr(tpsa, 'groups')
	if(missing(match.order)) {
		group.sizes <- table(tpsa$treat)
		ordering <- order(group.sizes, decreasing=FALSE)
		ordering <- ordering[c(2,1,3)]
		match.order <- names(group.sizes)[ordering]
	} else {
		ordering <- c(
			which(groups == match.order[1]),
			which(groups == match.order[2]),
			which(groups == match.order[3]))
	}

	ps1 <- getPS(tpsa, match.order[1], match.order[2])
	ps2 <- getPS(tpsa, match.order[2], match.order[3])
	ps3 <- getPS(tpsa, match.order[3], match.order[1])

	sd1 <- sd(ps1, na.rm=TRUE)
	sd2 <- sd(ps2, na.rm=TRUE)
	sd3 <- sd(ps3, na.rm=TRUE)
	
	if(!missing(exact)) {
		if(class(exact) %in% c('matrix','data.frame')) {
			exact <- apply(exact, 1, paste, collapse='.')
		}
	
		results <- data.frame()
		for(i in unique(exact)) {
			tmp <- tpsa[exact == i,]
			if(all(table(tmp$treat) > 1)) {
				results <- rbind(results, trimatch.apply2(
					tmp, caliper=caliper, nmatch=nmatch, match.order=match.order,
					sd1=sd1, sd2=sd2, sd3=sd3
				))
			} #Else there are not at least one unit from each group that exactly match
		}
	} else {
		results <- trimatch.apply2(tpsa, caliper, nmatch, match.order=match.order,
								   sd1=sd1, sd2=sd2, sd3=sd3)
	}
	
	names(results) <- c(match.order,
			paste0('D.', getModel(tpsa, match.order[1],match.order[2])),
			paste0('D.', getModel(tpsa, match.order[2],match.order[3])),
			paste0('D.', getModel(tpsa, match.order[3],match.order[1]))
	)
	
	results$Dtotal <- results$D.m1 + results$D.m2 + results$D.m3
	results <- results[order(results$Dtotal),]
	
	if(!is.null(method)) {
		results <- method(results, ...)
		row.names(results) <- 1:nrow(results)
	}
	
	unmatched <- tpsa[!(tpsa$id %in% c(results[,1], results[,2], results[,3])),]

	class(results) <- c('triangle.matches','data.frame')
	attr(results, 'triangle.psa') <- tpsa
	attr(results, 'match.order') <- match.order
	attr(results, 'unmatched') <- unmatched
	attr(results, getModel(tpsa, match.order[1],match.order[2])) <- c(match.order[1],match.order[2])
	attr(results, getModel(tpsa, match.order[2],match.order[3])) <- c(match.order[2],match.order[3])
	attr(results, getModel(tpsa, match.order[3],match.order[1])) <- c(match.order[3],match.order[1])
	
	return(results)
}

getModel <- function(tpsa, x, y) {
	test <- tpsa[which(tpsa$treat %in% c(x,y)),c('model1','model2','model3')]
	if(length(which(is.na(test$model1))) == 0) {
		return('m1')
	} else if(length(which(is.na(test$model2))) == 0) {
		return('m2')
	} else if(length(which(is.na(test$model3))) == 0) {
		return('m3')
	} else {
		stop('Could not find model.')
	}
}

getPS <- function(tpsa, x, y) {
	test <- tpsa[which(tpsa$treat %in% c(x,y)),c('model1','model2','model3')]
	if(length(which(is.na(test$model1))) == 0) {
		return(tpsa$ps1)
	} else if(length(which(is.na(test$model2))) == 0) {
		return(tpsa$ps2)
	} else if(length(which(is.na(test$model3))) == 0) {
		return(tpsa$ps3)
	} else {
		stop('Could not find model.')
	}
}

#' This method will return at least one treatment from groups one and two within
#' the caliper.
#' 
#' This method will attempt to return enough rows to use each treatment (the first
#' two groups in the matching order) at least once. Assuming treat1 is the first
#' group in the match order and treat2 the second, all duplicate treat1 rows
#' are removed. Next, all treat2 units not in present in after removing duplicate
#' treat1 units are identified. For each of those treat2 units, the matched
#' triplet with the smallest overall distances where treat2 is one of the mathched
#' units is retained.
#' 
#' @param tmatch initial results from \code{\link{trimatch}} that contains all
#'        possible matches within the specified caliper.
#' @param ... currently unused.
#' @export
maximumTreat <- function(tmatch, ...) {
	keep <- which(!duplicated(tmatch[,1]))
	t2.ids <- unique(tmatch[keep,2])
	retain <- tmatch[keep,]
	keep.t2 <- tmatch[-which(tmatch[,2] %in% t2.ids),]
	keep.t2 <- keep.t2[!duplicated(keep.t2[,2]),]
	retain <- rbind(retain, keep.t2)
	return(retain)
}

#' This method will use a M1-to-M2-to-1 matching.
#' 
#' In this method, \code{M2} corresponds to the number of times a treat1 unit can be
#' matched with a treat2 unit. The \code{M1} parameter corresponds to the number of
#' times a treat1 unit can be used in total.
#' 
#' @param tmatch initial results from \code{\link{trimatch}} that contains all
#'        possible matches within the specified caliper.
#' @param M1 a scaler indicating the number of unique subjects in group one to
#'        retain. This applies only to the first group in the matching order.
#' @param M2 a scaler indicating the number of unique matches to retain. This applies
#'        to the first two groups in the matching order.
#' @param ... currently unused.
#' @export
OneToN <- function(tmatch, M1=2, M2=1, ...) {
	results <- tmatch
	retain <- data.frame()
	for(i in 1:M2) {
		keep <- which(!duplicated(results[, 2]))
		if(length(keep) > 0) {
			retain <- rbind(retain, results[keep,])
			results <- results[-keep,]
		}
	}
	results <- retain
	row.names(results) <- 1:nrow(results)
	
	retain <- data.frame()
	for(i in 1:M1) {
		keep <- which(!duplicated(results[,1]))
		if(length(keep) > 0) {
			retain <- rbind(retain, results[keep,])
			results <- results[-keep,]
		}
	}
	return(retain)
}

#' Euclidean distance calculation.
#' 
#' This method uses a simple Euclidean distance caluclation for determining the
#' distances between two matches. That is, |ps1 - ps2|.
#' 
#' @param x vector of propensity scores.
#' @param grouping vector or factor identifying group membership.
#' @param id vector corresponding to unique identifer for each element in 
#'        \code{x} and \code{grouping}.
#' @param groups vector of length two indicating the unique groups to calculate
#'        the distance between. The first element will be the rows, the second columns.
#' @param caliper a scaler indicating the caliper to use for matching within
#'        each step.
#' @param nmatch number of smallest distances to retain.
#' @return a list of length equal to \code{x}. Each element of the list is a
#'        named numeric vector where the values correspond to the distance and the
#'        name to the \code{id}.
distance.euclid <- function(x, grouping, id, groups, caliper, nmatch=Inf) {
	sd <- sd(x, na.rm=TRUE)
	tmp <- data.frame(ps=x, group=grouping, id=id, stringsAsFactors=FALSE)
	tmp <- tmp[!is.na(tmp$ps),]
	tmp.1 <- tmp[tmp$group == groups[1],]
	tmp.2 <- tmp[tmp$group == groups[2],]
	d.list <- list()
	for(i in 1:nrow(tmp.1)) {
		dist <- abs(tmp.1[i,]$ps - tmp.2$ps)
		names(dist) <- tmp.2$id
		dist <- dist[order(dist)]
		dist <- dist / sd #Convert to standized units
		dist <- dist[1:min(length(dist), nmatch)]
		dist <- dist[dist < caliper]
		d.list[[i]] <- dist
	}
	names(d.list) <- as.character(tmp.1$id)
	return(d.list)
}
