#' Prints the results of \code{\link{mlpsa}} as a LaTeX table.
#' 
#' This function implements the \code{\link{xtable}} method for \code{\link{mlpsa}}.
#' 
#' @param x results of \code{\link{mlpsa}}
#' @param digits number of digits to print.
#' @param include.note include a table note indicating how many rows were removed
#'        due to insufficient data within a strata.
#' @param caption passed through to \code{\link{xtable}}.
#' @param label passed through to \code{\link{xtable}}.
#' @param align Not used. 
#' @param display passed through to \code{\link{xtable}}.
#' @param auto passed through to \code{\link{xtable}}.
#' @param ... other parameters passed to \code{\link{summary.mlpsa}}
#' @method xtable mlpsa
#' @import xtable
#' @export
xtable.mlpsa <- function(x, caption, label, align, digits=2, display=NULL,
						 auto=FALSE, include.note=TRUE, ...) {
	xdf <- suppressMessages(summary(x, ...))
	
	xdf$ci = NA
	xdf[!is.na(xdf$ci.min) & !is.na(xdf$ci.max),]$ci <- paste('(', 
			format(xdf[!is.na(xdf$ci.min) & !is.na(xdf$ci.max),]$ci.min, digits=digits), ', ', 
			format(xdf[!is.na(xdf$ci.min) & !is.na(xdf$ci.max),]$ci.max, digits=digits), ')', 
			sep='')
	xdf$ci.min <- NULL
	xdf$ci.max <- NULL
	xdf$Treat.n <- as.integer(xdf$Treat.n)
	xdf$Control.n <- as.integer(xdf$Control.n)
	xtab <- xtable(xdf, digits=digits, caption=caption, label=label, display=display, auto=auto,
			align=c('l','l','r','r','r@{\\extracolsep{.25cm}}','r','r','c'), ...)
	class(xtab) <- c('xmlpsa', class(xtab))
	attr(xtab, 'treat.name') <- x$x.label
	attr(xtab, 'control.name') <- x$y.label
	if(length(x$removed) > 0) {
		attr(xtab, 'note') <- paste('\\textit{Note:} Removed ', length(x$removed), ' (', 
				format(100*(length(x$removed) / x$overall.n), digits=1), 
				'\\%) rows due to insufficent data within a given strata.', sep='')
	}
	return(xtab)
}

#' Prints the results of \code{\link{mlpsa}} and \code{\link{xtable.mlpsa}}.
#' 
#' Print method for \code{\link{xtable.mlpsa}}.
#' 
#' @param x result of \code{\link{xtable.mlpsa}}
#' @param tabular.environment see \code{\link{print.xtable}}.
#' @param floating see \code{\link{print.xtable}}.
#' @param ... other parameters passed to \code{\link{print.xtable}}
#' @method print xmlpsa
#' @export
print.xmlpsa <- function(x, tabular.environment='longtable', floating=FALSE, ...) {
	xdf <- as.data.frame(x)
	hlineafter <- which(!is.na(xdf$level2)) - 1
	addtorow = list()
	addtorow$pos = list()
	addtorow$pos[[1]] = c(0)
	addtorow$command = c(paste(
		'\\hline & & \\multicolumn{2}{c}{', 
		attr(x, 'treat.name'), '} & \\multicolumn{2}{c}{', 
		attr(x, 'control.name'), '} & \\\\ \\cline{3-4} \\cline{5-6} ',
		'Level & Strata & Mean & n & Mean & n & Confidence Interval \\\\ \\hline', 
		'\\endfirsthead \\multicolumn{7}{l}{{...continued from previous page}}\\\\ \\hline ',
		' & & \\multicolumn{2}{c}{', 
		attr(x, 'treat.name'), '} & \\multicolumn{2}{c}{', 
		attr(x, 'control.name'), '} & \\\\ \\cline{3-4} \\cline{5-6} ',
		'Level & Strata & Mean & n & Mean & n & Confidence Interval \\\\ \\hline', 
		' \\endhead \\endfoot \\endlastfoot ',
		sep=''))
	if(!is.null(attr(x, 'note'))) {
		addtorow$pos[[2]] = c(nrow(x))
		addtorow$command = c(addtorow$command,
				paste('\\hline \\multicolumn{7}{l}{', attr(x, 'note'), '} \\\\', sep=''))
	}
	print.xtable(x, 
				tabular.environment=tabular.environment,
				floating=floating,
				include.rownames=FALSE, include.colnames=FALSE,
				add.to.row=addtorow,
				hline.after=hlineafter, ...)
}

