#' Convert \code{sim_reg_summary} object to \code{data.frame}
#'
#' @param x Object of class `sim_reg_summary`.
#' @param ... Non-used arguments.
#' @export
as.data.frame.sim_reg_summary <- function(x, ...) {
	top3 <- sort(x[["phi_marginal_term_freqs"]], decreasing=TRUE)[1:min(length(x[["phi_marginal_term_freqs"]]), 3)]
	data.frame(
		`P(gamma=1|y)`=x[["mean_posterior_gamma"]],
		`Top 3 phi terms`=paste(names(top3), ": ", round(digits=2, top3), sep="", collapse=", "),
		check.names=FALSE
	)
}

#' Print \code{sim_reg_summary} object
#'
#' @param x Object of class \code{sim_reg_summary}.
#' @param width Integer value for width of output
#' @param ontology Optional `ontology_index` class object allowing additional information about terms sampled in `phi` to be printed.
#' @param ... Non-used arguments.
#' @export
print.sim_reg_summary <- function(x, width=getOption("width"), ontology=NULL, ...) {
	numeric_params <- Filter(f=Negate(is.na), x=sim_reg_parameters[sim_reg_parameters_gamma & sim_reg_parameters_type == "numeric"])

	dashed <- paste0(rep("-", width), collapse="")
	cat(dashed, "\n")
	cat("P(gamma=1|y) = ", mean(x[["mean_posterior_gamma"]]), "\n")

	cat(dashed, "\n")
	if (x[["mean_posterior_gamma"]] > 0) {
		cat("Numeric parameters:\n")
		print(data.frame(
			Parameter=numeric_params,
			Mean=round(digits=2, x[["numeric_param_means_given_gamma1"]]),
			SD=round(digits=2, x[["numeric_param_sds_given_gamma1"]])
		), row.names=FALSE)

		cat(dashed, "\n")

		cat("Phi:\n")
		print(data.frame(
			t=names(x[["phi_marginal_term_freqs"]]),
			Name=if (!is.null(ontology)) substr(ontology$name[names(x[["phi_marginal_term_freqs"]])], 1, max(10, as.integer(width)-20)) else "-",
			P=round(x[["phi_marginal_term_freqs"]], digits=2),
			row.names=names(x[["phi_marginal_term_freqs"]])
		)[1:min(length(x[["phi_marginal_term_freqs"]]), 10),,drop=FALSE], row.names=FALSE)

		cat(dashed, "\n")
	}
}

#' Print \code{sim_reg_samples} object
#'
#' @param x Object of class \code{sim_reg_samples}.
#' @param width Integer value for width of output
#' @param ontology Optional \code{ontology_index} class object allowing additional information about terms sampled in \code{phi} to be printed.
#' @param ... Non-used arguments.
#' @export
print.sim_reg_samples <- function(x, width=getOption("width"), ontology=NULL, ...) {
	stopifnot(class(x) == "sim_reg_samples")
	print.sim_reg_summary(summary(x))
}

#' Get summary of \code{sim_reg_samples} object.
#'
#' @param object Object of class \code{sim_reg_samples}.
#' @param ... Non-used arguments
#' @export
#' @importFrom stats setNames
summary.sim_reg_samples <- function(object, ...) {
	numeric_params <- Filter(f=Negate(is.na), x=sim_reg_parameters[sim_reg_parameters_gamma & sim_reg_parameters_type == "numeric"])
	margs.tab <- sort(term_marginals(object$phi[object$gamma]), decreasing=TRUE)
	margs <- setNames(nm=names(margs.tab), as.numeric(margs.tab))

	structure(
		class="sim_reg_summary",
		list(
			mean_posterior_gamma=object[["mean_posterior_gamma"]],
			numeric_param_means_given_gamma1=sapply(object[numeric_params], mean),
			numeric_param_sds_given_gamma1=sapply(object[numeric_params], sd),
			phi_marginal_term_freqs=margs
		)
	)
}

#' Plot \code{sim_reg_summary} object
#'
#' @param x Object of class \code{sim_reg_summary}.
#' @template ontology
#' @param max_plot_terms Maximum number of terms (with non-zero sample frequency) to show in plot.
#' @param cex Text size for plot.
#' @param ... Other arguments to pass to \code{onto_plot}.
#' @export
#' @importFrom ontologyPlot onto_plot remove_links calibrate_sizes official_labels
#' @importFrom ontologyIndex get_ancestors
#' @importFrom grDevices hsv
#' @importFrom graphics lines
plot.sim_reg_summary <- function(x, ontology, max_plot_terms=10, cex=1, ...) {
	layout(matrix(c(1,2,2),ncol=3,nrow=1))

	lines.required <- length(x[["numeric_param_means_given_gamma1"]]) + 2

	plot(x=NULL, xlim=0:1, ylim=0:1, axes=FALSE, xlab="", ylab="")
	
	text(adj=c(0, 0.5), x=0.1, y=1-1/lines.required, labels=paste("P(gamma=1|y) = ", x[["mean_posterior_gamma"]], sep=""), cex=cex)
	lines(x=c(0.1, 0.9), y=rep(times=2, 1-2/lines.required))

	for (param.no in 1:length(x[["numeric_param_means_given_gamma1"]])) {
		text(adj=c(0, 0.5), x=0.1, y=1-(param.no+2)/lines.required, labels=paste(names(x[["numeric_param_means_given_gamma1"]])[param.no], " = ", round(digits=2, x[["numeric_param_means_given_gamma1"]][param.no]), sep=""), cex=cex)
	}

	subject.terms <- names(sort(decreasing=TRUE, x[["phi_marginal_term_freqs"]])[1:min(length(x[["phi_marginal_term_freqs"]]), max_plot_terms)])
	terms <- remove_links(ontology, get_ancestors(ontology, subject.terms))
	plot(onto_plot(
		ontology=ontology, 
		terms=terms, 
		fontsize=30,
		width=calibrate_sizes(ifelse(terms %in% subject.terms, x[["phi_marginal_term_freqs"]][terms], 0), 3, 1),
		label=ifelse(terms %in% subject.terms, paste(sep="\n", official_labels(ontology, terms), round(digits=2, x[["phi_marginal_term_freqs"]][terms])), official_labels(ontology, terms)),
		fillcolor=ifelse(terms %in% subject.terms, hsv(h=1/3, s=0.5), "light grey"),
		...
	))
}
