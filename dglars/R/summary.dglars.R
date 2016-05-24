summary.dglars <- function(object,k=c("BIC","AIC"),complexity=c("df","gdf"),digits = max(3, getOption("digits") - 3),...){
	if(is.numeric(k)){
		if(k<=0) stop("k must be greater than zero")
		knm <- "GoF"
	}
	else{
		knm <- match.arg(k)
		k <- ifelse(knm == "BIC",log(dim(object$X)[1]),2)
	}
	complexity <- match.arg(complexity)
	if(object$family=="poisson" & complexity == "gdf"){
		complexity <- "df"
		warning("'complexity' was set equal to 'df' because for Poisson regression the\n model complexity can be approximated by the number of nonzero coefficients")
	}
	tbl <- make_summary_table(object,k,complexity)
	names(tbl$table)[4] <- complexity
	names(tbl$table)[5] <- knm
	action <- object$action
	id <- which(action!="")
	n.tbl <- dim(tbl$table)[1]
	n.space <- length(id)
	id.tbl <- vector(length=n.tbl+n.space,mode="numeric")
	id.space <- id + seq(1,n.space)
	id.tbl[-id.space] <- seq(1:n.tbl)
	id.tbl[id.space] <- id
	tbl.format <- format(tbl$table[id.tbl,], digits = digits)
	tbl.format[id.space-1,1] <- ""
	tbl.format[id.space,-1] <- ""
	b.gof.names <- names(tbl$b.gof)[-1]
	if(is.null(b.gof.names)) b.gof.names <- "1"
	bestmodel.formula <- paste("y",paste(b.gof.names,collapse=" + "),sep=" ~ ")
	cat("\nCall:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print.data.frame(tbl.format, print.gap = 2, quote = FALSE, row.names = FALSE, ...)
	cat("\n==============================================\n")
	cat("\nBest model identified by",knm,"criterion ( k =",k,"and complexity =",complexity,"):\n\n", bestmodel.formula,"\n")
	cat("\nCoefficients:\n\n")
	print.default(format(tbl$b.gof, digits = digits), print.gap = 2, quote = FALSE,...)
	cat("\n\n",knm,": ",format(min(tbl$table[5]), digits = digits))
	cat("\n\n===\n\nAlgorithm", object$control$algorithm,"( method =",object$control$method,") with exit =",object$conv,"\n\n")
	invisible(tbl)
}
