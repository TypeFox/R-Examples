# Author: Robert J. Hijmans
# Date :  December 2009
# Version 0.1
# Licence GPL v3


evaluateROCR <- function(model, p, a, x) {
	if (! requireNamespace('ROCR')) {
		stop('ROCR package not found')
	}
	
	if (!missing(x)) {
		p <- predict(model, data.frame(extract(x, p)))
		a <- predict(model, data.frame(extract(x, a)))
	} else if (is.vector(p) & is.vector(a)) {
			# do nothing
	} else {
		p <- predict(model, data.frame(p))
		a <- predict(model, data.frame(a))
	}
	p <- stats::na.omit(p)
	a <- stats::na.omit(a)
	if (length(p) < 1) { stop('no valid presence (p) values') }
	if (length(a) < 1) { stop('no valid absence (a) values') }
	predictions = c(p, a)
	labels = c( rep(1, length(p)), rep(0, length(a)) )
	pred <- ROCR::prediction( predictions, labels)
	return(pred)
}




.auctest <- function(e) {
	w <- wilcox.test(e@presence, e@absence)
	pauc <- w$p.value
	w$auc <- as.vector(w$statistic) / (e@na * e@np)
}



evaluate <- function(p, a, model, x, tr, ...) {
	if (! missing(x) ) {
		p <- predict(model, data.frame(extract(x, p)), ...)
		a <- predict(model, data.frame(extract(x, a)), ...)
	} else if (is.vector(p) & is.vector(a)) {
			# do nothing
	} else {
		p <- predict(model, data.frame(p), ...)
		a <- predict(model, data.frame(a), ...)
	}
	p <- stats::na.omit(p)
	a <- stats::na.omit(a)
	np <- length(p)
	na <- length(a)
	if (na == 0 | np == 0) {
		stop('cannot evaluate a model without absence and presence data that are not NA')
	}

	if (missing(tr)) {
		if (length(p) > 1000) {
			tr <- as.vector(quantile(p, 0:1000/1000))
		} else {
			tr <- p
		}
		if (length(a) > 1000) {
			tr <- c(tr, as.vector(quantile(a, 0:1000/1000)))
		} else {
			tr <- c(tr, a)
		}
		tr <- sort(unique( round(tr, 8)))
		tr <- c( tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
	} else {
		tr <- sort(as.vector(tr))
	}
	
	N <- na + np

	xc <- new('ModelEvaluation')
	xc@presence = p
	xc@absence = a
		
	R <- sum(rank(c(p, a))[1:np]) - (np*(np+1)/2)
	xc@auc <- R / (as.numeric(na) * as.numeric(np))
	
	cr <- try( cor.test(c(p,a), c(rep(1, length(p)), rep(0, length(a))) ), silent=TRUE )
	if (class(cr) != 'try-error') {
		xc@cor <- cr$estimate
		xc@pcor <- cr$p.value
	}
	
	res <- matrix(ncol=4, nrow=length(tr))
	colnames(res) <- c('tp', 'fp', 'fn', 'tn')
	xc@t <- tr
	for (i in 1:length(tr)) {
		res[i,1] <- length(p[p>=tr[i]])  # a  true positives
		res[i,2] <- length(a[a>=tr[i]])  # b  false positives
		res[i,3] <- length(p[p<tr[i]])    # c  false negatives
		res[i,4] <- length(a[a<tr[i]])    # d  true negatives
	}
	xc@confusion = res
	a = res[,1]
	b = res[,2]
	c = res[,3]
	d = res[,4]
# after Fielding and Bell	
	xc@np <- as.integer(np)
	xc@na <- as.integer(na)
	xc@prevalence = (a + c) / N
	xc@ODP = (b + d) / N
	xc@CCR = (a + d) / N
	xc@TPR = a / (a + c)
	xc@TNR = d / (b + d)
	xc@FPR = b / (b + d)
	xc@FNR = c/(a + c)
	xc@PPP = a/(a + b)
	xc@NPP = d/(c + d)
	xc@MCR = (b + c)/N
	xc@OR = (a*d)/(c*b)

	prA = (a+d)/N
	prY = (a+b)/N * (a+c)/N
	prN = (c+d)/N * (b+d)/N
	prE = prY + prN
	xc@kappa = (prA - prE) / (1-prE)
	return(xc)
}




setMethod ('show' , 'ModelEvaluation', 
	function(object) {
		cat('class          :' , class(object), '\n')
		cat('n presences    :' , object@np, '\n')
		cat('n absences     :' , object@na, '\n')
		cat('AUC            :' , object@auc,'\n')
#		cat('p(AUC)      :' , object@pauc,'\n')
		cat('cor            :' , object@cor,'\n')
#		cat('p(cor)      :' , object@pcor,'\n')
#		cat('prevalence     :' , object@prevalence,'\n')
#		cat('overallDiagnosticPower :', object@ODP,'\n')
#		cat('correctClassificationRate :', object@CCR,'\n')
#		cat('sensitivity    :', object@TPR,'\n')
#		cat('specificity    :', object@TNR,'\n')
#		cat('falsePositiveRate :', object@FPR,'\n')
#		cat('falseNegativeRate :', object@FNR,'\n')
#		cat('PPP :', object@PPP,'\n')
#		cat('NPP :', object@NPP,'\n')
#		cat('misclassificationRate :', object@MCR,'\n')
#		cat('oddsRatio :', object@OR,'\n')
#		cat('kappa :', object@kappa,'\n')
		cat('max TPR+TNR at :', object@t[which.max(object@TPR + object@TNR)], '\n')
	}
)	


