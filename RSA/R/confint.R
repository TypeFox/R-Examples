#' @title Computes confidence intervals for RSA parameters, standard or bootstrapped
#' @description
#' Computes confidence intervals for RSA parameters, standard or bootstrapped (using a percentile bootstrap)
#'
#' @details
#' There are two ways of getting bootstrapped CIs and p-values in the RSA package If you use the option \code{se="boot"} in the \code{\link{RSA}} function, \code{lavaan} provides CIs and p-values based on the bootstrapped standard error (\emph{not} percentile bootstraps). If you use \code{confint(..., method="boot")}, in contrast, you get CIs and p-values based on percentile bootstrap.
#'
#' @method confint RSA
#' @export
#' @aliases confint

#' @param object An RSA object
#' @param parm Not used.
#' @param level The confidence level required.
#' @param digits Number of digits the output is rounded to; if NA, digits are unconstrained
#' @param model A string specifying the model; defaults to "full"
#' @param method "standard" returns the CI for the lavaan object as it was computed. "boot" computes new percentile bootstrapped CIs.
#' @param R If \code{method = "boot"}, R specifies the number of bootstrap samples
#' @param ... Additional parameters passed to the bootstrapLavaan function, e.g., \code{parallel="multicore", ncpus=2}.
#'
#' @seealso \code{\link{RSA}}
#'
#' @examples
#'
#' set.seed(0xBEEF)
#' n <- 300
#' err <- 2
#' x <- rnorm(n, 0, 5)
#' y <- rnorm(n, 0, 5)
#' df <- data.frame(x, y)
#' df <- within(df, {
#' 	diff <- x-y
#' 	absdiff <- abs(x-y)
#' 	SD <- (x-y)^2
#' 	z.sq <- SD + rnorm(n, 0, err)
#' })
#' 
#' r1 <- RSA(z.sq~x*y, df, models="SSQD")
#' (c1 <- confint(r1, model="SSQD"))
#'
#' # Dummy example with 10 bootstrap replications - better use >= 5000!
#' (c2 <- confint(r1, model="SSQD", method="boot", R=10))
#' \dontrun{
#' # multicore version
#' confint(r1, model="SSQD", R=5000, parallel="multicore", ncpus=2)
#' }


confint.RSA <- function(object, parm, level = 0.95, ..., model = "full", digits=3, method="standard", R = 5000) {	
	method <- match.arg(method, c("standard", "boot"))
	
	if (method == "standard") {
		p1 <- data.frame(parameterEstimates(object$models[[model]], level=level))
		p1 <- p1[p1$label != "", ]
		rownames(p1) <- p1$label
		p1 <- p1[, c("ci.lower", "ci.upper", "pvalue")]
		colnames(p1)[1:2] <- paste0(c((1-level)/2, 1-(1-level)/2)*100, "%")
		attr(p1, "type") <- paste0("CIs extracted from lavaan object.")
		if (!is.na(digits)) p1 <- round(p1, digits)
		return(p1)
	}
	
	if (method == "boot") {
		print(paste0("Drawing ", R, " bootstrap samples, please be patient ..."))

		# run parameterEstimates once to get the labels of the relevant indeces
		p1 <- data.frame(parameterEstimates(object$models[[model]]))
		
		b1 <- data.frame(bootstrapLavaan(object$models[[model]], FUN=function(y) {
				c1 <- coef(y, type="user")
				c1 <- c1[names(c1) %in% p1$label]
				return(c1)
			} , R = R, ...))
		
		CIs <- data.frame()
		for (i in 1:ncol(b1)) {
			x <- b1[, i]
			if (colnames(b1)[i] != "S") {
				# default: test against 0
				p <- sum(x<0)/length(x)	
			} else {
				# Parameter S: test against 1 (this is the null)
				p <- sum(x<1)/length(x)	
				message("Note: p value for parameter S indicates whether 1 (not zero!) is inside the CI.")
			}
			
			qu <- quantile(x, probs=c((1-level)/2, 1-(1-level)/2))
			res <- c(qu, pvalue=min(p, 1-p)*2)	# *2 to make p-values two-sided
			CIs <- rbind(CIs, res)
		}
	
		rownames(CIs) <- colnames(b1)
		colnames(CIs) <- c("LL", "UL", "pvalue")
		attr(CIs, "type") <- paste0("Bootstrapped CIs from ", R, " replications.\n", 
			ifelse("S" %in% colnames(b1), "Note: p value for parameter S indicates whether 1 (not zero!) is inside the CI.", "")
		)
		if (!is.na(digits)) CIs <- round(CIs, digits)
		return(CIs)
	}	
}


