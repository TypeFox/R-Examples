#' Summary of mpqtl object
#' 
#' Prints a summary of the detected QTL
#' @S3method summary mpqtl
#' @method summary mpqtl
#' @param object Object of class \code{mpqtl}
#' @param ... Additional arguments
#' @return Table with rows for each QTL detected:
#' Column 1 is the chromosome where the QTL was detected
#' Column 2 is the position where the QTL was detected on the chromosome
#' Columns 3 and 4 are the flanking markers for the QTL
#' Columns 5, 6, 7 and 8 are the effect estimates for the founders
#' Column 9 is the Wald test statistic for the overall test at that position
#' Column 10 is the p-value for the test statistic
#' @seealso \code{\link[mpMap]{mpIM}}, \code{\link[mpMap]{plot.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, thr=1, responsename="pheno")
#' summary(mpq.dat)

summary.mpqtl <- function(object, ...)
{
 	n.founders <- nrow(object$founders)
	map <- object$map
	f3 <- substr(rownames(object$founders), 1, 3)
	fmap <- attr(object$prob, "map")
	qtlres <- object$QTLresults$qtl
	nqtl <- attr(qtlres, "nqtl")
	if (nqtl==0) cat("No QTL found.\n")
	else 
	{
		chr <- rep(names(qtlres), unlist(lapply(qtlres, function(x) return(nrow(x)))))
		cm <- round(unlist(lapply(qtlres, function(x) return(x[,1]))),2)
		effect <- round(matrix(unlist(lapply(qtlres, function(x) return(as.vector(t(x[, 1+1:n.founders]))))), nrow=n.founders, ncol=length(cm)), 3)
		se <- round(matrix(unlist(lapply(qtlres, function(x) return(as.vector(t(x[, 1+n.founders+1:n.founders]))))), nrow=n.founders, ncol=length(cm)), 3)

		fmrkl <- vector(length=nqtl)
		fmrkr <- vector(length=nqtl)
		wald <- vector(length=nqtl)
		pval <- vector(length=nqtl)
		degf <- vector(length=nqtl)
		index <- vector(length=nqtl)
		for (ii in 1:nqtl)
		{
			mrkli <- which.max(map[[chr[ii]]]*(map[[chr[ii]]]<=cm[ii]))
			if (length(map[[chr[ii]]]) > 1) {
			if (mrkli==length(map[[chr[ii]]])) mrkli <- mrkli-1
			fmrkl[ii] <- names(map[[chr[ii]]])[mrkli]
			fmrkr[ii] <- names(map[[chr[ii]]])[mrkli+1]
			} else fmrkl[ii] <- fmrkr[ii] <- names(map[[chr[ii]]])[mrkli]
			wald[ii] <- round(object$QTLresults$wald[[chr[ii]]][which.min(abs(cm[ii]-fmap[[chr[ii]]]))], 2)
			pval[ii] <- signif(object$QTLresults$pval[[chr[ii]]][which.min(abs(cm[ii]-fmap[[chr[ii]]]))], 3)
			degf[ii] <- signif(object$QTLresults$degf[[chr[ii]]][which.min(abs(cm[ii]-fmap[[chr[ii]]]))], 3)
		}
		eff3 <- paste("Effect_",f3,sep="")
		se3 <- paste("SE_",f3,sep="")
	if (n.founders==4)
		table <- data.frame("Chr"=chr, "Pos"=cm, "LeftMrk"=fmrkl, "RightMrk"=fmrkr, effect[1,], se[1,], effect[2,], se[2,], effect[3,], se[3,], effect[4,], se[4,], "Wald"=wald, "df"=degf, "pvalue"=pval)
	else if (n.founders==8)
		table <- data.frame("Chr"=chr, "Pos"=cm, "LeftMrk"=fmrkl, "RightMrk"=fmrkr, effect[1,], se[1,], effect[2,], se[2,], effect[3,], se[3,], effect[4,], se[4,], effect[5,], se[5,], effect[6,], se[6,], effect[7,], se[7,], effect[8,], se[8,], "Wald"=wald, "df"=degf, "pvalue"=pval)

		names(table)[seq(5, 4+2*n.founders, 2)] <- eff3
		names(table)[seq(6, 4+2*n.founders, 2)] <- se3

		table
	}
}

