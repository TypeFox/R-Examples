#' @method print fSRM
#' @export

print.fSRM <-
function(x, digits=3, ..., var.onesided=TRUE) {
	
	# Print model summary, also show package version
	if(!exists("meta") || is.null(meta)) meta <- packageDescription("fSRM")
	cat(sprintf("fSRM version %s", meta$Version))
	cat("\n================================\n\n")
	## The model for 4 members must have 31 free parameters and 47 df!
	cat(paste("SRM with roles (Roles: ", paste(x$roles, collapse=", "), "); DVs = ", paste0(x$var.id, collapse=", "), "\n", sep=""))
	if (!is.null(x$group)) cat(paste0("Split by group variable `", x$group,"`\n"))
	
	
	cat("\nModel summary:\n----------------\n")
	show(x$fit)
	cat("\nModel Fit:\n----------------\n")
	FIT <- fitmeasures(x$fit)
	cat(paste("Chi2 (df=", FIT["df"], ") = ", round(FIT["chisq"], digits), ", p = ", round(FIT["pvalue"], digits), "\n", sep=""))
	cat(paste("CFI = ", round(FIT["cfi"], digits), "\n", sep=""))
	cat(paste("TLI / NNFI = ", round(FIT["tli"], digits), "\n", sep=""))
	cat(paste("RMSEA = ", round(FIT["rmsea"], digits), " [", round(FIT["rmsea.ci.lower"], digits), ";", round(FIT["rmsea.ci.upper"], digits), "]", "; Test of close fit: p(data | true value == .05) = ", round(FIT["rmsea.pvalue"], digits), "\n", sep=""))
	
	
	if (is.null(x$group)) {
		print.singlegroup(x, group=1, digits=digits, var.onesided=var.onesided)
	} else {
		# print stats for each group
		for (g in 1:length(x$groupnames)) {
			cat("\n\n#####################################\n")
			cat(paste0("Statistics for group", g, "\n"))
			cat("#####################################\n")
			print.singlegroup(x, group=g, digits=digits, var.onesided=var.onesided)
		}
		
		if (x$diff == TRUE & x$means == TRUE) {
			cat("\n\n#####################################\n")
			cat(paste0("Difference of means between groups (", x$groupnames[1], "-", x$groupnames[2], ")\n"))
			cat("#####################################\n")
		
			MD <- getCor(x, label=".meanDiff.", group=0)[, -1]			
			colnames(MD)[1] <- "component"
			colnames(MD)[2] <- "diff"
			
			pval <- MD$p.value
			MD2 <- cbind(MD[, 1:5], sig=p2star(pval), MD[, 6:7])
			MD2$p.value <- p(pval, digits)
			print(MD2, row.names=FALSE)
		}
		if (x$diff == TRUE) {
			cat("\n\n#####################################\n")
			cat(paste0("Difference of variances between groups (", x$groupnames[1], "-", x$groupnames[2], ")\n"))
			cat("#####################################\n")
		
			VD <- getCor(x, label=".varDiff.", group=0)[, -1]
			colnames(VD)[1] <- "component"
			colnames(VD)[2] <- "diff"
			
			pval <- VD$p.value
			VD2 <- cbind(VD[, 1:5], sig=p2star(pval), VD[, 6:7])
			VD2$p.value <- p(pval, digits)
			print(VD2, row.names=FALSE)
		}
	}
}






print.singlegroup <-
function(x, group=1, digits=3, conf.level=0.95, var.onesided=TRUE) {
	
	if (var.onesided == TRUE) {
		conf.level <- 1-(1-conf.level)*2
	}
	
	eff <- as.data.frame(parameterEstimates(x$fit, level=conf.level))
	eff$f <- paste(eff$lhs, eff$op, eff$rhs)
	
	# SS = standardized solution: get correlation for that
	SS <- getCor(x, ops=c("~~", "~"), group=group)

	cat("\n\nVariance decomposition:\n----------------\n")
	T <- varComp(x, group=group, conf.level=conf.level)
	if (var.onesided == TRUE) {
		T$p.value <- T$p.value/2
	}
	pval <- T$p.value
	T[, -1] <- round(T[, -1], digits)
	T2 <- cbind(T[, 1:5], sig=p2star(pval), T[, 6:7])
	T2$p.value <- p(pval, digits)
	print(T2)
	if (var.onesided == TRUE) {
		cat("\n(p-values are for one-sided tests for variances; confidence level for CIs is", round(conf.level*100, 2), "%)\n\n")
	}
	
	
	cat("\n\nRelative variance decomposition:\n----------------\n")
	print(round(percTable(x, group=group)$stand * 100))
	

	if (!x$drop %in% c("actor", "partner", "GR")) {
		cat("\n\nGeneralized reciprocity (actor-partner covariances):\n----------------\n")
		GR <- getGR(x, group=group)
		pval <- GR$p.value
		GR2 <- cbind(GR[, 1:5], sig=p2star(pval), GR[, 6:8])
		GR2$p.value <- p(pval, digits)
		GR2$r <- f2(GR2$r, digits=digits, skipZero=TRUE)
		print(GR2)
	}
	
	#cat("\n\nDyadic reciprocity (relationship covariances): Mean r =", round(meanNA(GR$COR), digits),"(out of bounds estimates set to NA)\n----------------\n")
	DR <- getDR(x, group=group)
	#cat("\n\nDyadic reciprocity (relationship covariances): Mean r =", round(meanNA(DR$r), digits),"(out of bounds estimates set to NA)\n----------------\n")
	cat("\n\nDyadic reciprocity (relationship covariances):\n----------------\n")
	DR2 <- cbind(DR[, 1:5], sig=p2star(DR$p.value), DR[, 6:8])
	DR2$p.value <- p(DR2$p.value, digits)
	DR2$r <- f2(DR2$r, digits=digits, skipZero=TRUE)
	print(DR2)
	
	if (length(x$IGSIM) > 0) {
		cat("\n\nIntragenerational similarity:\n----------------\n")
		igsim <- SS[grepl("IGSIM", SS$label), ][, -2]
		pval <- igsim$p.value
		igsim2 <- cbind(igsim[, 1:5], sig=p2star(pval), igsim[, 6:8])
		igsim2$p.value <- p(igsim2$p.value, digits)
		igsim2$r <- f2(igsim2$r, digits=digits, skipZero=TRUE)
		print(igsim2)
	}
	
	# TODO: Include self-ratings
	# if (x$self == TRUE) {
# 		AS <- data.frame()
# 		for (t in x$roles) {
#
# 			if (x$selfmode == "cor") {F <- paste(style$actor, ".", t, " ~~ ", style$self, ".", t, sep="")}
# 			if (x$selfmode == "kq") {F <- paste(style$self, ".", t, " ~ ", style$actor, ".", t, sep="")}
#
# 			AS0 <- SS[SS$f == F, ]
# 			AS0$comment <- ""
#
# 			# get Variance of components --> if that is < min.p, correlation is not reliable!
# 			SD1 <- SS[SS$f == paste(style$partner, ".", t, " ~~ ", style$partner, ".", t, sep=""), ]
# 			SD2 <- SS[SS$f == paste(style$self, ".", t, " ~~ ", style$self, ".", t, sep=""), ]
# 			if (is.na(SD1$pvalue)) SD1$pvalue <- 1
# 			if (is.na(SD2$pvalue)) SD2$pvalue <- 1
#
# 			if (SD1["pvalue"] > x$min.p | SD2["pvalue"] > x$min.p) {
# 				AS0$COR <- NA_real_
# 				AS0$comment <- paste("One of the variance components has p <", x$min.p)
# 			}
#
# 			if (AS0$pvalue > x$min.p) {
# 				AS0$COR <- NA_real_
# 				AS0$comment <- paste("Covariance estimate has p <", x$min.p)
# 			}
#
# 			AS <- rbind(AS, AS0)
# 		}
# 		#cat("\n\nAssumed similarity: Mean r =", round(meanNA(AS$COR), digits),"(out of bound estimates set to zero)\n----------------\n")
# 		cat("\n\nAssumed similarity:\n----------------\n")
# 		print(AS,row.names=TRUE)
#
# 		SO <- data.frame()
# 		for (t in x$roles) {
# 			if (x$selfmode == "cor") {F <- paste(style$partner, ".", t, " ~~ ", style$self, ".", t, sep="")}
# 			if (x$selfmode == "kq") {F <- paste(style$self, ".", t, " ~ ", style$partner, ".", t, sep="")}
#
# 			SO0 <- SS[SS$f == F, ]
# 			SO0$comment <- ""
#
# 			# get Variance of components --> if that is < min.p, correlation is not reliable!
# 			SD1 <- SS[SS$f == paste(style$partner, ".", t, " ~~ ", style$partner, ".", t, sep=""), ]
# 			SD2 <- SS[SS$f == paste(style$self, ".", t, " ~~ ", style$self, ".", t, sep=""), ]
# 			if (is.na(SD1$pvalue)) SD1$pvalue <- 1
# 			if (is.na(SD2$pvalue)) SD2$pvalue <- 1
#
# 			if (SD1["pvalue"] > x$min.p | SD2["pvalue"] > x$min.p) {
# 				SO0$COR <- NA_real_
# 				SO0$comment <- paste("One of the variance components has p <", x$min.p)
# 			}
#
# 			if (SO0$pvalue > x$min.p) {
# 				SO0$COR <- NA_real_
# 				SO0$comment <- paste("Covariance estimate has p <", x$min.p)
# 			}
#
# 			SO <- rbind(SO, SO0)
# 		}
# 		SO$COR <- as.numeric(SO$COR)
# 		#cat("\n\nSelf-Other agreement: Mean r =", round(meanNA(SO$COR), digits),"(out of bound estimates set to NA)\n----------------\n")
# 		cat("\n\nSelf-Other agreement:\n----------------\n")
# 		print(SO, row.names=TRUE)
# 	}
	
	if (x$means == TRUE) {
		if (x$pairwise==TRUE) {
			cat("\n\nMean structure: Indices starting with 'C.' are pairwise comparisons between roles\n----------------\n")
		} else {
			cat("\n\nMean structure\n----------------\n")
		}
		if (is.null(x$group)) {
			MS <- eff[grepl(".means.", eff$label, fixed=TRUE), c(1, 5:10)]
		} else {
			MS <- eff[grepl(paste0(".means", x$groupnames[group], "."), eff$label, fixed=TRUE), c(1, 6:11)]
		}
		colnames(MS) <- c("factor", "estimate", "se", "z", "p.value", "ci.lower", "ci.upper")		
		rownames(MS) <- NULL

		pval <- MS$p.value
		MS[, -1] <- round(MS[, -1], digits)
		MS2 <- cbind(MS[, 1:5], sig=p2star(pval), MS[, 6:7])		
		MS2$p.value <- p(pval, digits)
		
		print(MS2)
	}
	
}
