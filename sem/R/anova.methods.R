# last modified 2012-08-01 by J. Fox

anova.objectiveML <- function(object, model.2, robust=FALSE, ...){
	anovaAdjchisq <- function(adjobj0, adjobj1){
		# this subfunction oringally by Jarrett Byrnes
		#from http://www.statmodel.com/chidiff.shtml
		# Satorra-bentler adjusted chi sq
		#switching to get order right
		sbs.nested <- adjobj0
		sbs.full <- adjobj1
		t0 <- sbs.nested$chisq
		tr0 <- sbs.nested$chisq.scaled
		t1 <- sbs.full$chisq
		tr1 <- sbs.full$chisq.scaled		
		c0 <- sbs.nested$c
		c1 <- sbs.full$c		
		d0 <- sbs.nested$df
		d1 <- sbs.full$df		
		cd <- (d0 * c0 - d1*c1)/(d0 - d1)
		trd <- abs((t0 - t1)/cd) 		
		df <- abs(d0 - d1)		
		table <- data.frame(c(d0, d1), c(tr0, tr1), c(NA, df), c(NA, trd),
				c(NA, pchisq(trd, df, lower.tail=FALSE)))		
		return(table)
	}
	dev.1 <- deviance(object)
	df.1 <- df.residual(object)
	dev.2 <- deviance(model.2)
	df.2 <- df.residual(model.2)
	name.1 <- deparse(substitute(object))
	name.2 <- deparse(substitute(model.2))
	df <- abs(df.1 - df.2)
	if (df == 0) stop("the models have the same Df")
	if (object$N != model.2$N)
		stop("the models are fit to different numbers of observations")
	if ((nrow(object$S) != nrow(model.2$S)) || !all.equal(object$S, model.2$S))
		stop("the models are fit to different moment matrices")
	if(!robust){
		chisq <- abs(dev.1 - dev.2)
		table <- data.frame(c(df.1, df.2), c(dev.1, dev.2), c(NA, df), c(NA, chisq),
				c(NA, pchisq(chisq, df, lower.tail=FALSE)))
	}
	else{
		cat("Adjusted Using Satorra-Bentler Correction\n");
		table <- anovaAdjchisq(object$adj.obj, model.2$adj.obj)
	}
	names(table) <- c("Model Df", "Model Chisq", "Df", "LR Chisq", "Pr(>Chisq)")
	rownames(table) <- c(name.1, name.2)
	structure(table, heading = c("LR Test for Difference Between Models", ""),
			class = c("anova", "data.frame"))
}

anova.objectiveFIML <- function(object, model.2, ...){
    logLik.1 <- logLik(object)
    df.1 <- df.residual(object)
    logLik.2 <- logLik(model.2)
    df.2 <- df.residual(model.2)
    name.1 <- deparse(substitute(object))
    name.2 <- deparse(substitute(model.2))
    df <- abs(df.1 - df.2)
    if (df == 0) stop("the models have the same Df")
    if (object$N != model.2$N)
        stop("the models are fit to different numbers of observations")
    if ((nrow(object$S) != nrow(model.2$S)) || !all.equal(object$S, model.2$S))
        stop("the models are fit to different data sets")
    chisq <- 2*(abs(logLik.1 - logLik.2))
    table <- data.frame(c(df.1, df.2), c(logLik.1, logLik.2), c(NA, df), c(NA, chisq),
                        c(NA, pchisq(chisq, df, lower.tail=FALSE)))
    names(table) <- c("Model Df", "Model Log-Likelihood", "Df", "LR Chisq", "Pr(>Chisq)")
    rownames(table) <- c(name.1, name.2)
    structure(table, heading = c("LR Test for Difference Between Models", ""),
              class = c("anova", "data.frame"))
}

