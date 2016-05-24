

#' A Simple Binomial Test for Type I Error in mvQCA
#'
#' A Binomial test for multi-value qualitative comparative analysis (mvQCA),
#' designed to calculate the probability of a false positive given the number of
#' hypotheses implicitly tested and the number of confirming cases.
#' @param freq.y The frequency with which the dependent variable occurs in the
#' sample (the number of 1s divided by the total number of cases).
#' @param configs A list of configurations and the number of cases in which each
#' configuration occurs.
#' @param total.configs The total number of configurations used in the original
#' mvQCA analysis. This will generally equal the number of lines in the
#' truth table used for Boolean minimization.
#' @param adj.method The method used to calculate adjusted p-values (see 
#' \code{\link{p.adjust}} for details).
#' @return An object containing the results of the Binomial test.
#' @keywords Binomial mvQCA multiple inference p-value adjust
#' @export
#' @examples
#' test <- mvQCAbinTest(freq.y=0.7, configs=list(aB=5, bCD=3, Ce=2),
#'    total.configs=20)
#' summary(test)


mvQCAbinTest <- function(freq.y, configs, total.configs, adj.method="holm"){
	if(typeof(configs)!="list"){stop("configs must be a list of the form list(name1=number.of.cases1, name2=number.of.cases2, ...)")}
	if(is.null(names(configs))){stop("configs must be a list of the form list(name1=number.of.cases1, name2=number.of.cases2, ...)")}
	if(total.configs < length(configs)){stop("Value of total.configs cannot be less than the length of configs.")}
	if(length(freq.y)>1 | freq.y<0 | freq.y>1){stop("y must be a single value between 0 and 1.")}
	if(sum(configs > (total.configs - (total.configs*freq.y)))>0){warning("Number of confirming cases not consistent with total number of configurations and frequency of outcome.", call.=FALSE)}

	out.mat <- matrix(, nrow = length(configs), ncol = 3)
	out.mat[,1] <- as.vector(unlist(configs))
	out.mat[,2] <- (1-freq.y)^(out.mat[,1])
	out.mat[,3] <- p.adjust(c(out.mat[,2], rep(1,(total.configs-length(configs)))), method=adj.method)[1:length(configs)]
	rownames(out.mat) <- names(configs)
	colnames(out.mat) <- c("Number of cases", "p-raw", "p-adj")
	return.obj <- list(call = match.call(), total.configurations = total.configs, config.names = names(configs), p.adj.method = adj.method, result = out.mat)
	class(return.obj) <- "mvQCAbt"
	return.obj
}
