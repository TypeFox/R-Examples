#' This function will perform phase II of the multilevel propensity score analysis.
#' 
#' TODO: Need more details
#' 
#' The ci.adjust provides a Bonferroni-Sidak adjusted confidence intervals based
#' on the number of levels/clusters.
#' 
#' @param response vector containing the response values
#' @param treatment vector containing the treatment conditions
#' @param strata vector containing the strata for each response
#' @param level2 vector containing the level 2 specifications
#' @param minN the minimum number of subjects per strata for that strata to be 
#'        included in the analysis.
#' @param reverse reverse the order of treatment and control for the difference
#'        calculation.
#' @param ci.level the confidence level to use for confidence intervals. Defaults
#'        to a 95\% confidence level.
#' @return a mlpsa class
#' @seealso \code{\link{mlpsa.ctree}} \code{\link{mlpsa.logistic}}
#' @export
#' @examples
#' \dontrun{
#' require(multilevelPSA)
#' require(party)
#' data(pisana)
#' data(pisa.colnames)
#' data(pisa.psa.cols)
#' mlctree = mlpsa.ctree(pisana[,c('CNT','PUBPRIV',pisa.psa.cols)], formula=PUBPRIV ~ ., level2='CNT')
#' student.party = getStrata(mlctree, pisana, level2='CNT')
#' student.party$mathscore = apply(student.party[,paste0('PV', 1:5, 'MATH')], 1, sum) / 5
#' results.psa.math = mlpsa(response=student.party$mathscore, 
#'        treatment=student.party$PUBPRIV, 
#'        strata=student.party$strata, 
#'        level2=student.party$CNT, minN=5)
#' results.psa.math
#' summary(results.psa.math)
#' }
mlpsa <- function(response, treatment=NULL, strata=NULL, level2=NULL, 
				  minN=5, reverse=FALSE, ci.level=0.05) {
	stopifnot(length(response) == length(treatment)  & 
		length(treatment) == length(strata) & length(strata) == length(level2))
	
	if(reverse) {
		x.pos <- 5
		y.pos <- 4
	} else {
		x.pos <- 4
		y.pos <- 5
	}
	
	multilevelPSA <- list()
	
	xscale <- .1
	yscale <- .1
	
	thedata <- data.frame(response=response, treatment=treatment, strata=strata, level2=level2)
	thedata$strata2 <- paste(thedata$level2, thedata$strata, sep='.')
	t <- as.data.frame(table(thedata$strata2, thedata$treatment, useNA='ifany'))
	#Remove groupings where there is not at least one record in each group (i.e. public and private)
	multilevelPSA$removed <- which(thedata$strata2 %in% t[which(t$Freq < minN), 'Var1'])
	thedata <- thedata[-multilevelPSA$removed,]
	message(paste('Removed ', (length(response) - nrow(thedata)), 
				  ' (', prettyNum(100 * (length(response) - nrow(thedata)) / 
				  length(response), digits=3), '%) rows due to stratum size being less than ', 
				  minN, sep=''))
	thedata$strata2 <- as.factor(as.character(thedata$strata2))
	thedata$level2 <- as.factor(as.character(thedata$level2))
	
	#Summary statistics by each stratum
	d <- describeBy(thedata$response, list(thedata$treatment, thedata$strata2), 
					mat=TRUE, skew=FALSE)
	d <- d[,c('group1', 'group2', 'n', 'mean', 'se')]
	names(d) <- c('treatment', 'strata2', 'n', 'Mean', 'se')
	d <- cbind(cast(d, strata2 ~ treatment, value='n'), cast(d, strata2 ~ treatment, value='Mean')[,2:3], cast(d, strata2 ~ treatment, value='se')[,2:3])
	names(d)[2] <- paste(names(d)[2], 'n', sep='.')
	names(d)[3] <- paste(names(d)[3], 'n', sep='.')
	names(d)[6] <- paste(names(d)[6], 'se', sep='.')
	names(d)[7] <- paste(names(d)[7], 'se', sep='.')
	d$n <- d[,y.pos-2] + d[,x.pos-2]
	d$Diff <- d[,y.pos] - d[,x.pos]
	mapping <- thedata[!duplicated(thedata$strata2),c('level2', 'strata2')]
	d <- merge(d, mapping, by='strata2', all.x=TRUE)
	
	#Weighted means by level 2
	diff.wtd <- data.frame(level2=character(), n=integer(), diffwtd=numeric(), 
						   mnx=numeric(), mny=numeric(), mnxy=numeric(), 
						   ci.min=numeric(), ci.max=numeric(), df=numeric())
	nlevels <- length(unique(d$level2))
	for(i in unique(d$level2)) {
		tmp <- d[which(d$level2==i),]
		n <- sum(tmp$n)
		wtss <- tmp$n/n
		mny <- sum(tmp[,y.pos] * wtss)
		mnx <- sum(tmp[,x.pos] * wtss)
		yn <- sum(tmp[,y.pos-2])
		xn <- sum(tmp[,x.pos-2])
		mnxy <- (mnx + mny)/2
		diffwtd <- sum(tmp$Diff * tmp$n) / n
		
		#Calculate confidence interval
		tmp <- thedata[which(thedata$level2 == i),]
		n <- length(tmp$response)
		nstrat <- dim(table(factor(tmp$strata2)))
		ncontrol <- as.data.frame(table(factor(tmp$strata2), tmp$treatment))[1:nstrat, 3]
		ntreat <- as.data.frame(table(factor(tmp$strata2), tmp$treatment))[(nstrat + 1):(2 * nstrat), 3]
		o <- order(tmp$treatment)
		ord.strata <- factor(tmp$strata2[o])
		nc <- table(tmp$treatment)[1]
		nt <- table(tmp$treatment)[2]
		ord.response <- tmp$response[o]
		var.0 <- tapply(ord.response[1:nc], ord.strata[1:nc], var)
		ni.0 <- table(ord.strata[1:nc])
		frac.0 <- var.0/ncontrol
		ncp1 <- nc + 1
		ncpnt <- nc + nt
		var.1 <- tapply(ord.response[ncp1:ncpnt], ord.strata[ncp1:ncpnt], var)
		ni.1 <- table(ord.strata[ncp1:ncpnt])
		frac.1 <- var.1/ntreat
		se.wtd <- ((sum(frac.0) + sum(frac.1))^0.5)/nstrat
		ci.diff <- diffwtd
		df <- length(tmp$response) - 2 * length(unique(tmp$strata2))
		ci.min <- ci.diff - qt((1 - (ci.level/2)), df) * se.wtd
		ci.max <- ci.diff + qt((1 - (ci.level/2)), df) * se.wtd
		ci.level.adjust <- (1-(1-ci.level)^(1/nlevels)) / 2
		ci.min.adjust <- ci.diff - qt(1-ci.level.adjust, df) * se.wtd
		ci.max.adjust <- ci.diff + qt(1-ci.level.adjust, df) * se.wtd
		
		diff.wtd <- rbind(diff.wtd, data.frame(level2=i, n=n, diffwtd=diffwtd, 
							mnx=mnx, mny=mny, mnxy=mnxy,
							xn=xn, yn=yn,
							ci.min=ci.min, ci.max=ci.max, 
							df=df, se.wtd=se.wtd,
							ci.min.adjust=ci.min.adjust, ci.max.adjust=ci.max.adjust))
	}
	names(diff.wtd)[which(names(diff.wtd)=='mnx')] <- names(d)[x.pos]
	names(diff.wtd)[which(names(diff.wtd)=='mny')] <- names(d)[y.pos]
	names(diff.wtd)[which(names(diff.wtd)=='xn')] <- names(d)[x.pos-2]
	names(diff.wtd)[which(names(diff.wtd)=='yn')] <- names(d)[y.pos-2]
	
	multilevelPSA$x.label <- names(d)[x.pos]
	multilevelPSA$y.label <- names(d)[y.pos]
	multilevelPSA$overall.n <- sum(d$n)
	multilevelPSA$overall.nx <- sum(d[,x.pos-2])
	multilevelPSA$overall.ny <- sum(d[,y.pos-2])
	multilevelPSA$overall.wtss <- d$n / multilevelPSA$overall.n
	multilevelPSA$overall.mnx <- sum(d[,x.pos] * multilevelPSA$overall.wtss)
	multilevelPSA$overall.mny <- sum(d[,y.pos] * multilevelPSA$overall.wtss)
	multilevelPSA$overall.mnxy <- (multilevelPSA$overall.mnx + multilevelPSA$overall.mny) / 2
	multilevelPSA$overall.wtd <- sum(d$Diff * d$n) / multilevelPSA$overall.n
	multilevelPSA$overall.se.wtd <- sum(diff.wtd$se.wtd * diff.wtd$n) / multilevelPSA$overall.n
	multilevelPSA$approx.t <- multilevelPSA$overall.wtd / multilevelPSA$overall.se.wtd
	
	#Calculate confidence interval (borrowed from circ.psa in PSAgraphics package)
	n <- length(thedata$response)
	nstrat <- dim(table(thedata$strata2))
	ncontrol <- as.data.frame(table(thedata$strata2, thedata$treatment))[1:nstrat, 3]
	ntreat <- as.data.frame(table(thedata$strata2, thedata$treatment))[(nstrat + 1):(2 * nstrat), 3]
	o <- order(thedata$treatment)
	ord.strata <- thedata$strata2[o]
	nc <- table(thedata$treatment)[1]
	nt <- table(thedata$treatment)[2]
	ord.response <- thedata$response[o]
	var.0 <- tapply(ord.response[1:nc], ord.strata[1:nc], var)
	ni.0 <- table(ord.strata[1:nc])
	frac.0 <- var.0/ncontrol
	ncp1 <- nc + 1
	ncpnt <- nc + nt
	var.1 <- tapply(ord.response[ncp1:ncpnt], ord.strata[ncp1:ncpnt], var)
	ni.1 <- table(ord.strata[ncp1:ncpnt])
	frac.1 <- var.1/ntreat
	se.wtd <- ((sum(frac.0) + sum(frac.1))^0.5)/nstrat
	ci.diff <- multilevelPSA$overall.wtd
	df <- length(thedata$response) - 2 * length(unique(thedata$level2))
	multilevelPSA$overall.ci <- c(ci.diff - qt(0.975, df) * se.wtd, ci.diff + qt(0.975, df) * se.wtd)
	
	#Unweighted means
	d2 <- describeBy(thedata$response, list(thedata$level2, thedata$treatment), mat=TRUE)
	d2 <- d2[,c('group1', 'group2', 'mean')]
	names(d2) <- c('level2', 'treatment', 'Mean')
	d2 <- cast(d2, level2 ~ treatment, value='Mean')
	
	multilevelPSA$plot.range <- c(min(d[,4:5]) - xscale * (max(d[,4:5]) - min(d[,4:5])),
			max(d[,4:5]) + yscale * (max(d[,4:5]) - min(d[,4:5])))
	#multilevelPSA$projection.intercept <- 2 * (min(d[,4:5]) - yscale * (max(d[,4:5]) - min(d[,4:5])))
	multilevelPSA$projection.intercept <- 2 * min(d[,4:5]) + (.1 + yscale) * (max(d[,4:5]) - min(d[,4:5]))
	
	diff.wtd$xmark <- (multilevelPSA$projection.intercept - diff.wtd$diffwtd) / 2
	diff.wtd$ymark <- diff.wtd$xmark + diff.wtd$diffwtd
	
	multilevelPSA$level1.summary <- d
	multilevelPSA$unweighted.summary <- as.data.frame(d2)
	multilevelPSA$level2.summary <- diff.wtd
	
	class(multilevelPSA) <- 'mlpsa'
	
	return(multilevelPSA)
}

