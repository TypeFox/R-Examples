#' Create summary of \code{BeviMed} classed-object
#'
#' @param object Object of class \code{BeviMed}.
#' @param ... Other arguments passed to \code{\link{estimate_confidence_interval}}.
#' @return Object of class \code{BeviMed_summary}.
#' @export
summary.BeviMed <- function(object, ...) {
	vt <- object[["variant_table"]]
	y <- object[["y"]]
	variant_counts <- lapply(setNames(nm=c(F,T)), function(y_is) as.integer(table(factor(vt$variant[y[vt$case]==y_is], levels=1:object[["k"]]))))

	structure(list(
		Z=apply(object[["Z"]][,object[["k"]]*(length(object[["temperatures"]])-1)+1:object[["k"]],drop=FALSE], 2, mean),
		phi=mean(exp(object[["log_phi"]])),
		omega=mean(1-1/(1+exp(object[["logit_omega"]]))),
		confidence_interval=estimate_confidence_interval(
			temperatures=object[["temperatures"]],
			y_log_lik_t_equals_1_traces=object[["y_log_lik_t_equals_1"]],
			...
		),
		ML_v=power_posteriors_ML_sum(object[["y_log_lik_t_equals_1"]], object[["temperatures"]]),
		ML_n=`log P_n`(object[["y"]]),
		samples=nrow(object[["y_log_lik_t_equals_1"]]),
		estimate_omega=object[["estimate_omega"]],
		estimate_phi=object[["estimate_phi"]],
		phi_acceptance_rate=apply(object[["log_phi"]], 2, function(log_phis) mean(log_phis[-length(log_phis)] != log_phis[-1])), 
		omega_acceptance_rate=apply(object[["logit_omega"]], 2, function(logit_omegas) mean(logit_omegas[-length(logit_omegas)] != logit_omegas[-1])), 
		n=length(object[["y"]]),
		k=object[["k"]],
		variant_counts=variant_counts,
		temperatures=object[["temperatures"]]
	), class="BeviMed_summary")
}

#' Print readable summary of \code{BeviMed_summary} object.
#'
#' @param x \code{BeviMed_summary} object.
#' @param print_Z Logical value indicating whether to print list of marginal probabilities of \code{Z_j = 1} for all variants \code{j}.
#' @param width Width of printing.
#' @param ... Not-used arguments 
#' @return Prints a summary
#' @export
print.BeviMed_summary <- function(x, print_Z=FALSE, width=getOption("width"), ...) {
	stopifnot(class(x) == "BeviMed_summary")
	dashed <- paste0(rep("-", width), collapse="")
	cat(dashed, "\n")
	cat("Log Bayes factor between variant level model and null model is ", round(x[["ML_v"]]-x[["ML_n"]], digits=2), "\n", sep="")
	cat("\n")
	cat("A confidence interval for the estimate is:\n")
	print(round(digits=2, x[["confidence_interval"]] - x[["ML_n"]]))

	cat(dashed, "\n")
	if (x[["estimate_omega"]]) {
		cat("Estimate of omega: ", round(digits=2, x[["omega"]]), "\n", sep="")
		cat("\tAcceptance rate in sequential chains: \n\t", paste0(collapse=" : ", round(digits=2, x[["omega_acceptance_rate"]])), "\n", sep="")
		if (x[["estimate_phi"]]) {
			cat("Estimate of phi: ", round(digits=2, x[["phi"]]), "\n", sep="")
			cat("\tAcceptance rate in sequential chains: \n\t", paste0(collapse=" : ", round(digits=2, x[["phi_acceptance_rate"]])), "\n", sep="")
		}
		cat(dashed, "\n")
	}

	#cat("Estimates by MCMC samples based on ", x[["samples"]], " samples in ", length(x[["temperatures"]]), " tempered chains", "\n", sep="")
	#cat(dashed, "\n")

	#only show this for the highest temperature...
	if (print_Z) {
		cat("Estimated probabilities of pathogenicity of individual variants\n")
		cat("(conditional on model v)\n\n")

		print(row.names=FALSE, data.frame(
			check.names=FALSE,
			stringsAsFactors=FALSE,
			Variant=if (is.null(names(x[["Z"]]))) 1:length(x[["Z"]]) else names(x[["Z"]]),
			Controls=x[["variant_counts"]][["FALSE"]],
			Cases=x[["variant_counts"]][["TRUE"]],
			`P(Z_j=1|y,V)`=round(digits=2, x[["Z"]]),
			`Bar Chart`=sapply(1:length(x[["Z"]]), function(j) paste0("[", paste0(collapse="", rep("=", as.integer(x[["Z"]][j]*20))), paste0(collapse="", rep(" ", 20-as.integer(x[["Z"]][j]*20))), "]"))))

	}
}

#' Print readable summary of \code{BeviMed} object.
#'
#' @param x \code{BeviMed_summary} object.
#' @param ... Not-used arguments 
#' @return Prints a summary
#' @export
print.BeviMed <- function(x, ...) {
	stopifnot(class(x) == "BeviMed")
	print.BeviMed_summary(summary(x), ...)
}


