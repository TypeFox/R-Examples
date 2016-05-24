#'Visually comparing shrinkage and LS estimates
#'
#'This function creates a plot (using \code{qplot()}) where the shrinkage
#'estimate appears on the horizontal axis and the LS estimate appears on the
#'vertical axis.
#'
#'
#'@param eb a matrix of random effects
#'@param ols a matrix of the OLS estimates found using \code{random_ls_coef}
#'@param identify the percentage of points to identify as unusual,
#'\code{FALSE} if you do not want the points identified.
#'@param silent logical: should the list of data frames used to make the plots
#' be supressed.
#'@param ... other arguments to be passed to \code{qplot()}
#'@author Adam Loy \email{loyad01@@gmail.com}
#'@export
#'@examples
#'
#'wages.fm1 <- lmer(lnw ~ exper + (exper | id), data = wages)
#'wages.sepLM <- adjust_lmList(lnw ~ exper | id, data = wages)
#'rancoef.eb <- coef(wages.fm1)$id
#'rancoef.ols <- coef(wages.sepLM)
#'compare_eb_ls(eb = rancoef.eb, ols = rancoef.ols, identify = 0.01)
#'
compare_eb_ls <- function(eb, ols, identify = FALSE, silent = TRUE, ...){
  unusual <- ids <- NULL # Make codetools happy
	ret <- NULL
	for(i in 1:dim(ols)[2]){
	p <- qplot(x = eb[,i], y = ols[,i], geom = "point", main = colnames(eb)[i], 
	xlab = "EB resid", ylab = "LS resid", ...) + 
			geom_abline(intercept = 0, slope = 1, linetype = I(2)) +
			geom_smooth(method = "lm", se = FALSE) 
		if(identify != FALSE){
			temp_eb <- eb[,i]
			attr(temp_eb, "names") <- rownames(eb)
			extreme <- identify_resid(eb = temp_eb, ols = ols[,i], identify = identify)
			dat <- cbind(eb = temp_eb, ols = ols[,i])
			dat <- as.data.frame(dat)
			dat$ids <- names(temp_eb)
			extreme <- merge(x = dat, y = extreme)
			# extreme <- cbind(dat, extreme)
			# extreme <- extreme[order(abs(extreme$residual), decreasing = TRUE), ]
			# n <- round(dim(extreme)[1] * identify)
			
			p <- p + geom_text(data = subset(extreme, unusual == TRUE), aes(x = eb, y = ols, label = ids, hjust=.5, vjust=1.5), size=I(3))
			#ret <- list(ret, subset(extreme, unusual == TRUE))
			ret <- list(ret, extreme)
		}
	oask <- grDevices::devAskNewPage(TRUE)
        on.exit(grDevices::devAskNewPage(oask))
	print(p)
	}
	
	if(silent == FALSE){
		ret <- list(ret[[1]][[2]], ret[[2]])
		return(ret)
	}
}
