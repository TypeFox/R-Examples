#' Replace term ID names of vector/array with full term names
#'
#' @template ontology
#' @param x Vector/array named by term IDs
#' @export
full_names <- function(ontology, x) {
	if (is.array(x)) 
		dimnames(x) <- lapply(dimnames(x), function(i) ontology$name[i])
	else
		names(x) <- ontology$name[names(x)]
	x
}

#' Get marginal frequencies of single terms in phi
#'
#' @template phi
#' @return Numeric vector of proportion of samples each term was included in phi.
#' @export
term_marginals <- function(phi) {
	x <- sort(decreasing=TRUE, table(unlist(phi)))
	x/length(phi)
}

#' Get marginal frequencies of single terms in phi
#'
#' @template term_descendancy_matrix
#' @param phi_vector_trace Delta trace output from \code{\link{sim_reg}} function, where the rows are the values of phi at each iteration
#' @return Numeric vector of proportion of samples each term was included in phi.
term_marginals_from_phi_vec <- function(term_descendancy_matrix, phi_vector_trace) 
	with(data=as.numeric(sort(decreasing=TRUE, table(phi_vector_trace[as_row_leaves(term_descendancy_matrix, phi_vector_trace)]))), expr=x/sum(x))
 
#' Get matrix of presence of term pairs
#'
#' @template phi
#' @param terms_to_use Specify the terms to appear in the plot.
#' @param symmetric Logical value determining whether output is a symmetric matrix (\code{TRUE}) or lower triangular one (\code{FALSE}).
#' @param n Maximum number of terms to include in the matrix.
#' @return Matrix of term-pair frequencies (where by cell i,j contains the frequency of inclusion in phi for term pair i,j).
#' @export
#' @importFrom stats setNames
#' @importFrom utils combn 
term_pair_marginals <- function(
	phi,
	terms_to_use=NULL,
	symmetric=FALSE,
	n=10
) {
	all.terms <- unique(unlist(phi))
	n <- min(n, length(all.terms))
	selected <- 
		if (is.null(terms_to_use))
			names(sort(decreasing=TRUE, table(unlist(phi)))[1:n])
		else
			terms_to_use

	selected.full.trace <- sapply(phi, function(x) setNames(selected %in% x, selected)) 
	selected.trace <- selected.full.trace[,apply(selected.full.trace, 2, function(y) sum(y) >= 2)]

	term.pairs <- t(combn(selected, 2))
	heatmap.mat <- matrix(0, n, n, dimnames=rep(list(selected), 2))
	heatmap.mat[term.pairs] <- apply(
		term.pairs,
		1,
		function(pair) sum(apply(
			selected.trace[pair,],
			2,
			all
		))/length(phi)
	)

	heatmap.mat <- heatmap.mat + t(heatmap.mat)

	if (!symmetric) heatmap.mat[lower.tri(heatmap.mat, diag=TRUE)] <- NA

	heatmap.mat
}

#' Create `ontology_plot' of phi
#'
#' Create plot of marginal frequency of individual terms in context of related terms in the ontology, ready for plotting to device or exporting to dot file.
#'
#' @template ontology
#' @template phi
#' @param max_terms Specify the maximum number of terms to appear in the plot.
#' @param min_frequency Threshold frequency for including terms in plot.
#' @param colour_gradient Logical indicating whether to colour terms in the plot according to their marginal frequencies (blue being the least frequent, yellow the most).
#' @param size_gradient Logical indicating whether to colour terms in the plot according to their marginal frequencies.
#' @param custom_labels Character vector of custom labels for terms (named by corresponding term IDs).
#' @param show_proportion Logical indicating whether to append the `inclusion in phi' rate to the labels in the terms.
#' @param fillcolor Vector of colours (named by HPO term IDs) for corresponding terms in plot.
#' @param ... Additional parameters to be passed to \code{onto_plot}.
#' @return Plots graph.
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom ontologyPlot remove_links simple_labels onto_plot
#' @importFrom ontologyIndex get_ancestors
phi_plot <- function(
	ontology, 
	phi,
	max_terms=10,
	min_frequency=0,
	colour_gradient=FALSE,
	size_gradient=TRUE,
	custom_labels=c(),
	show_proportion=TRUE,
	fillcolor=NULL,
	...
) {
	marginal_freqs <- local({ x <- term_marginals(phi); y <- x[1:min(length(x), max_terms)]; y[y > min_frequency] })
	plot_term_freqs <- local({
		missing <- setdiff(remove_links(ontology, get_ancestors(ontology, names(marginal_freqs))), names(marginal_freqs))
		setNames(
			c(marginal_freqs, rep(0, length(missing))),
			c(names(marginal_freqs), missing)
		)
	})

	onto_plot(
		ontology,
		terms=names(plot_term_freqs), 
		label=
			if (show_proportion) paste(ifelse(!(names(plot_term_freqs) %in% names(custom_labels)), simple_labels(ontology, names(plot_term_freqs)), custom_labels[names(plot_term_freqs)]), paste(round(plot_term_freqs, 2), sep=""), sep="\n")
			else ifelse(!(names(plot_term_freqs) %in% names(custom_labels)), simple_labels(ontology, names(plot_term_freqs)), custom_labels[names(plot_term_freqs)])
		,
		width="if"(
			size_gradient,
			sqrt(plot_term_freqs)*(3/max(sqrt(plot_term_freqs))),
			rep(1, length(plot_term_freqs))
		),
		fillcolor=if (is.null(fillcolor)) "if"(
			colour_gradient,
			colorRampPalette(c("#0099FF", "green3", "Yellow"))(20)[local({ x <- cut(plot_term_freqs, breaks=seq(from=min(plot_term_freqs), to=max(plot_term_freqs), by=diff(range(plot_term_freqs))/20), include.lowest=TRUE, labels=FALSE); x[which(plot_term_freqs == max(plot_term_freqs))] <- 20; x })],
			rep("cyan", length(plot_term_freqs))
		) else fillcolor,
		...
	)
}
 
#' Plot marginal distribution of term-pair inclusion in phi
#'
#' Plots the marginal frequency of terms pairs' inclusion in the phi parameter over the course of an application of \code{\link{sim_reg}} as heatmap.
#'
#' @param x Numeric matrix.
#' @param ... Other parameters to pass to image.
#' @return Nothing - side effect: plots graph.
#' @export
#' @importFrom plotrix color.legend
#' @importFrom graphics image axis
#' @importFrom grDevices rgb 
term_pair_marginals_plot <- function(x, ...) {
	image(x, col=rgb(c(rep(0, 50), seq(from=0, to=1, by=0.02)), seq(from=0, to=1, by=0.01), seq(from=0.5, to=0, by=-0.005)), axes=FALSE, ...)
	axis(1, at=0:(ncol(x)-1)/(ncol(x)-1), las=2, labels=colnames(x))
	axis(2, las=2, labels=rownames(x), at=0:(nrow(x)-1)/(nrow(x)-1)) 
	color.legend(align="rb", gradient="y", rect.col=rgb(c(rep(0, 50), seq(from=0, to=1, by=0.02)), seq(from=0, to=1, by=0.01), seq(from=0.5, to=0, by=-0.005)), xl=1.02, xr=1.06, yb=0, yt=1, legend=as.character(signif(digits=2, range(x, na.rm=TRUE))))
}

#' Plot SimReg output
#'
#' Plot the output of the call to the \code{sim_reg} function. Output contains plots of marginal probabilities of terms and pairs of terms being present phi and histograms giving the distributions of variables. Note, a large device is required for a successful call.
#'
#' @param x Output of call to \code{sim_reg}.
#' @template ontology
#' @param ... Non-used arguments.
#' @return Plots
#' @export
#' @importFrom graphics curve hist layout plot.new text plot par legend
#' @importFrom stats dnorm pbeta qnorm 
plot.sim_reg_samples <- function(x, ontology, ...) {
	sim_reg_out <- x

	pbeta2 <- function(x, mean, alpha.plus.beta) pbeta(q=x, shape1=mean*alpha.plus.beta, shape2=(1-mean)*alpha.plus.beta)

	equal.w.cols <- c(2, 9)
	its <- length(sim_reg_out$gamma)

	layout(do.call(what=cbind, lapply(1:length(equal.w.cols), function(col.ind) sort(rep(Reduce(init=1, f="+", x=c(0, equal.w.cols)[1:col.ind]) + 0:(equal.w.cols[col.ind]-1), as.integer(200/equal.w.cols[col.ind])+1)[1:200]))))

	par(mar=c(5, 4, 4, 2))

	variables <- c("gamma", "alpha_star", "alpha", "log_beta", "logit_mean_f", "log_alpha_plus_beta_f", "logit_mean_g", "log_alpha_plus_beta_g", "phi")
	numeric.vars <- Filter(x=variables, f=function(x) class(sim_reg_out[[x]]) == "numeric")
	non.char <- setdiff(variables, "phi")
	no.gamma <- setdiff(variables, "gamma")

	plot.new()

	text(cex=3, "Mean parameter values", x=0, y=1, adj=0)
	par(mar=c(5, 4, 4, 2))
	for (i in 1:length(non.char)) {
		text(cex=2, adj=0, x=0, y=1-i/(2*length(non.char)), paste(non.char[i], ": ", round(digits=2, mean(sim_reg_out[[non.char[i]]])), sep=""))
	}

	for (i in 1:length(no.gamma)) {
		text(cex=2, adj=0, x=0, y=0.5-i/(2*length(no.gamma)), paste(no.gamma[i], " acceptance: ", round(digits=2, mean(sim_reg_out[[paste(no.gamma[i], "_accept", sep="")]])), sep=""))
	}
	par(mar=c(5, 4, 4, 2))

	if (sum(sim_reg_out$gamma) > 0) {
		plot(phi_plot(ontology, sim_reg_out$phi[sim_reg_out$gamma], font.size=100, max_terms=10, colour_gradient=FALSE), main=expression(paste("Marginal Posterior of Inclusion of Terms in ", phi, sep="")), cex.main=3)
	} else {
		plot.new()
	}

	par(mar=c(5, 4, 4, 2))

	alpha_star.range <- range(sim_reg_out$alpha_star, qnorm(p=c(0.05, 1-0.05), mean=sim_reg_out$priors$alpha_star_mean, sd=sim_reg_out$priors$alpha_star_sd))
	hist(main=NULL, x=sim_reg_out$alpha_star[sim_reg_out$gamma], ylab="Density", xlab="alpha star", freq=FALSE, breaks=seq(from=alpha_star.range[1], to=alpha_star.range[2], by=diff(alpha_star.range)/100), xlim=alpha_star.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim_reg_out$priors$alpha_star_mean, sd=sim_reg_out$priors$alpha_star_sd), col="blue")

	alpha.range <- range(sim_reg_out$alpha, qnorm(p=c(0.05, 1-0.05), mean=sim_reg_out$priors$alpha_mean, sd=sim_reg_out$priors$alpha_sd))
	hist(main=NULL, x=sim_reg_out$alpha[sim_reg_out$gamma], ylab="Density", xlab="alpha ", freq=FALSE, breaks=seq(from=alpha.range[1], to=alpha.range[2], by=diff(alpha.range)/100), xlim=alpha.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim_reg_out$priors$alpha_mean, sd=sim_reg_out$priors$alpha_sd), col="blue")

	log_beta.range <- range(sim_reg_out$log_beta, qnorm(p=c(0.05, 1-0.05), mean=sim_reg_out$priors$log_beta_mean, sd=sim_reg_out$priors$log_beta_sd))
	hist(main=NULL, x=sim_reg_out$log_beta[sim_reg_out$gamma], ylab="Density", xlab="log beta", freq=FALSE, breaks=seq(from=log_beta.range[1], to=log_beta.range[2], by=diff(log_beta.range)/100), xlim=log_beta.range)
	curve(add=TRUE, expr=dnorm(x, mean=sim_reg_out$priors$log_beta_mean, sd=sim_reg_out$priors$log_beta_sd), col="blue")

	hist(main=NULL, x=1-1/(1+exp(sim_reg_out$logit_mean_f[sim_reg_out$gamma])), xlim=c(0, 1), breaks=seq(from=0, to=1, by=0.01), freq=FALSE, xlab="mean f")
	curve(add=TRUE, col="blue", expr=dnorm(log(x)-log(1-x), mean=sim_reg_out$priors$logit_mean_f_mean, sd=sim_reg_out$priors$logit_mean_f_sd)/x/(1-x))

	log_alpha_plus_beta_f.range <- range(sim_reg_out$log_alpha_plus_beta_f[sim_reg_out$gamma], qnorm(p=c(0.05, 1-0.05), mean=sim_reg_out$priors$log_alpha_plus_beta_f_mean, sd=sim_reg_out$priors$log_alpha_plus_beta_f_sd))
	hist(main=NULL, x=sim_reg_out$log_alpha_plus_beta_f[sim_reg_out$gamma], breaks=seq(from=log_alpha_plus_beta_f.range[1], to=log_alpha_plus_beta_f.range[2], by=diff(log_alpha_plus_beta_f.range)/100), freq=FALSE, xlab="log_alpha_plus_beta_f")
	curve(add=TRUE, expr=dnorm(x, mean=sim_reg_out$priors$log_alpha_plus_beta_f_mean, sd=sim_reg_out$priors$log_alpha_plus_beta_f_sd), col="blue")

	plot(main="f transformation sample", x=NULL, xlim=c(0, 1), ylim=c(0, 1), ylab="f(similarity)", xlab="similarity")
	for (i in sample(which(sim_reg_out$gamma), size=min(sum(sim_reg_out$gamma), 100))) curve(add=T, col=rgb(0,0,1,0.3), expr=pbeta2(x, mean=1-1/(1+exp(sim_reg_out$logit_mean_f[i])), alpha.plus.beta=exp(sim_reg_out$log_alpha_plus_beta_f[i])))

	hist(main=NULL, x=1-1/(1+exp(sim_reg_out$logit_mean_g[sim_reg_out$gamma])), xlim=c(0, 1), breaks=seq(from=0, to=1, by=0.01), freq=FALSE, xlab="mean g")
	curve(add=TRUE, col="blue", expr=dnorm(log(x)-log(1-x))/x/(1-x))

	log_alpha_plus_beta_g.range <- range(sim_reg_out$log_alpha_plus_beta_g[sim_reg_out$gamma], qnorm(p=c(0.05, 1-0.05), mean=sim_reg_out$priors$log_alpha_plus_beta_g_mean, sd=sim_reg_out$priors$log_alpha_plus_beta_g_sd))
	hist(main=NULL, x=sim_reg_out$log_alpha_plus_beta_g[sim_reg_out$gamma], breaks=seq(from=log_alpha_plus_beta_g.range[1], to=log_alpha_plus_beta_g.range[2], by=diff(log_alpha_plus_beta_g.range)/100), freq=FALSE, xlab="log_alpha_plus_beta_g")
	curve(add=TRUE, expr=dnorm(x, mean=sim_reg_out$priors$log_alpha_plus_beta_g_mean, sd=sim_reg_out$priors$log_alpha_plus_beta_g_sd), col="blue")

	plot(main="g transformation sample", x=NULL, xlim=c(0, 1), ylim=c(0, 1), ylab="g(similarity)", xlab="similarity")
	for (i in sample(which(sim_reg_out$gamma), size=min(sum(sim_reg_out$gamma), 100))) curve(add=T, col=rgb(0,0,1,0.3), expr=pbeta2(x, mean=1-1/(1+exp(sim_reg_out$logit_mean_g[i])), alpha.plus.beta=exp(sim_reg_out$log_alpha_plus_beta_g[i])))

}

#' Create grid of SimReg output plots
#'
#' Create a PDF containing a set of summary plots for the output of the \code{sim_reg} application. Output contains plots of marginal probabilities of terms and pairs of terms being present phi and histograms giving the distributions of variables.
#'
#' @template ontology
#' @param file_name File to write plots to 
#' @param sim_reg_out Output of call to \code{sim_reg}
#' @return Plots graph to file
#' @export
#' @importFrom grDevices dev.off pdf
sim_reg_summary <- function(ontology, file_name, sim_reg_out) {
	pdf(file_name, width=25, height=35)
	plot(sim_reg_out, ontology=ontology)
	dev.off()
}
