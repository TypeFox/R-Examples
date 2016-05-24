#' @title Compare a full list of RSA models
#'
#' @description
#' Compare several fit indexes of all models computed from the RSA function
#'
#' @details
#' No details so far.
#'
#' @export
#' @param x An RSA object
#' @param verbose Should the summary be printed?
#' @param plot Should the comparison be plotted (using the \code{\link{modeltree}} function)?
#' @param digits Digits of the output
#' @param ... Additional parameters passed to the \code{\link{modeltree}} function


compare <- function(x, verbose=TRUE, plot=FALSE, digits=3, ...) {
	
	if (table(sapply(x$models, is.null))["FALSE"] <= 1) {
		stop("You need more than one models for comparison!")
	}
	
	with(x$models, {
	
	res <- data.frame()
	
	if (!is.null(full)) {
		free.max <- getFreeParameters(full)
	
		if (verbose==TRUE) {
			cat("-------------------------------------------------------------------------\n")	
			cat("Standard polynomial models:\n")
			cat("-------------------------------------------------------------------------\n")
		}
	
		res1 <- cModels(mL=list(cubic=cubic, full=full, IA=IA, additive=additive, diff=diff, null=null), set="directed", free.max=free.max)
		if (verbose==TRUE & !is.null(res1)) {
			cat("Testing directed difference models: Interaction, additive main effects, difference model :\n")
			cat("-------------------------------------------------------------------------\n")
			print(round(res1[, 1:13], digits))
		}
			
		res2 <- cModels(list(cubic=cubic, full=full, SRRR=SRRR, SRSQD=SRSQD, SSQD=SSQD, SQD=SQD, null=null), set="flat_sq", free.max)
		if (verbose==TRUE & !is.null(res2)) {
			cat("\n\nTesting 'flat ridge' discrepancy models against SRRR and full polynomial model:\n")
			cat("-------------------------------------------------------------------------\n")
			print(round(res2[, 1:13], digits))
		}
	
		res3 <- cModels(list(cubic=cubic, full=full, SRRR=SRRR, SRR=SRR, RR=RR, SQD=SQD, null=null), set="RR", free.max)
		if (verbose==TRUE & !is.null(res3)) {
			cat("\n\nTesting 'rising ridge' against full polynomial model:\n")
			cat("-------------------------------------------------------------------------\n")
			print(round(res3[, 1:13], digits))
		}
		
		## compute additional comparisons
		res4 <- cModels(list(full=full, SRR=SRR, SSQD=SSQD), set="SRR_SSQD", free.max)
		if (verbose==TRUE & !is.null(res4)) {
			cat("\n\nTesting transition from SRR to SSQD model (i.e., removing the mean level effect from SRR):\n")
			cat("-------------------------------------------------------------------------\n")
			print(round(res4[, 1:13], digits))
		}
		
		## single variable models (only x + x2, or y + y2)
		res5 <- cModels(list(full=full, onlyx2=onlyx2, onlyx=onlyx), set="onlyx", free.max)
		if (verbose==TRUE & !is.null(res5)) {
			cat("\n\nSingle variable models (only x + x^2):\n")
			cat("-------------------------------------------------------------------------\n")
			print(round(res5[, 1:13], digits))
		}
		res6 <- cModels(list(full=full, onlyy2=onlyy2, onlyy=onlyy), set="onlyy", free.max)
		if (verbose==TRUE & !is.null(res6)) {
			cat("\n\nSingle variable models (only y + y^2):\n")
			cat("-------------------------------------------------------------------------\n")
			print(round(res6[, 1:13], digits))
		}
		
		
		reslist <- list(res1, res2, res3, res4, res5, res6)
		res <- plyr::rbind.fill(reslist[!sapply(reslist, is.null)])
	}
		
		
	
	aL3 <- anovaList(list(absunc=absunc, absdiff=absdiff))
	if (aL3$n.mods > 1) {
		if (verbose==TRUE) {
			cat("\n\n-------------------------------------------------------------------------\n")	
			cat("Piecewise regression: absolute difference vs. unrestricted difference model\n")
			cat("-------------------------------------------------------------------------\n")
		}
		free.max2 <- getFreeParameters(absunc)
		a3 <- cbind(aL3$ANOVA, plyr::ldply(aL3$models, function(X) {
			F <- fitmeasures(X)
			R <- inspect(X, "r2")
			names(R) <- "R2"
			n <- lavaan::nobs(X)
			k <- free.max2 - F["df"]
			R2.p <- pf(((n-k-1)*R)/(k*(1-R)), k, n-k-1, lower.tail=FALSE)
			names(R2.p) <- "R2.p"
			return(c(F[c("cfi", "tli", "rmsea", "srmr")], R, R2.p))
		}))
		a3 <- a3[, !grepl(".id", colnames(a3))]
		a3$k <- free.max2 - a3$Df
		a3$R2.adj <- 1 - ((1-a3$R2))*((lavaan::nobs(absunc)-1)/(lavaan::nobs(absunc)-a3$k-1))
		a3$delta.R2 <- c(NA, a3$R2[1:(nrow(a3)-1)] - a3$R2[2:(nrow(a3))])
		if (verbose==TRUE) print(round(a3, digits))
		a3$model <- rownames(a3)
		a3$set <- "abs"
		res <- rbind(res, a3)
	}
	
	class(res) <- c("data.frame", "cRSA")
	
	if (plot==TRUE) {
		modeltree(res, ...)
	}
	
	invisible(res)
	})
}




#' @title Compare two specific RSA models
#'
#' @description
#' Compare several fit indexes of two models computed from the RSA function
#'
#' @details
#' You must take care yourself that the compared models are nested! There is no automatic check, so you could, in principle, compare non-nested models. This is valid for AIC, BIC, CFI, and R2 indices, but *not* for the chi2-LR test!
#'
#' @export
#' @param x An RSA object
#' @param m1 Name of first model
#' @param m2 Name of second model
#' @param digits Digits of the output
#' @param verbose Should the summary be printed?
compare2 <- function(x, m1="", m2="full", digits=3, verbose=TRUE) {

	if (is.null(x$models[[m1]]) | is.null(x$models[[m2]])) {
		stop("You need two models for comparison! At least one of the models has not been fit.")
	}
	

	free.max <- getFreeParameters(x$models[[m1]])
	mL <- list(M1=x$models[[m1]], M2=x$models[[m2]])
	names(mL) <- c(m1, m2)
	res <- cModels(mL, set="two_models", free.max)
	if (verbose==TRUE & !is.null(res)) {
		print(round(res[, 1:13], digits))
	}

	invisible(res)
}