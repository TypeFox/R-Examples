#' Subset \code{sim_reg_samples} object
#'
#' @param samples \code{sim_reg_samples} classed object.
#' @param at Integer vector of indices of iterations to extract from samples.
#' @export
subset_sim_reg_samples <- function(samples, at) {
	structure(class="sim_reg_samples", c(
		samples[setdiff(names(samples), sim_reg_all_traces)],
		lapply(FUN=function(trace) { if (length(trace) == 0) { trace } else { if (is.matrix(trace)) trace[at,] else trace[at] } }, X=samples[sim_reg_all_traces])
	))
}

#' Prune between-term similarity matrix
#'
#' Prune similarity of more specific terms to 0
#'
#' @template ontology
#' @template term_sim_mat
#' @return Numeric matrix of term-term similarities
#' @export
prune_sim_mat <- function(ontology, term_sim_mat) {
	result <- term_sim_mat
	for (term in rownames(result)) 
		result[term, ] <- ifelse(colnames(result) %in% ontology$ancestors[[term]], result[term,], 0)
	result
}

#' Get leaf matrix
#'
#' Procure logical matrix from character matrix of ontological terms, indicating whether each element is a leaf in the context of the row. Typically used to map the sampled vectors of the ontological parameter phi from the \code{\link{sim_reg}} procedure.
#'
#' @param term_descendancy_matrix Logical term_descendancy_matrix, dimensions symmetrically labelled by terms, and where by a cell value of TRUE indicates that the row is the ancestor of the column term (in the sense of the DAG structure of the HPO)
#' @param terms_matrix The character matrix of HPO terms
#' @return Logical matrix the same dimensions as the terms_matrix which indicates whethere each element is a leaf or not
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' as_row_leaves(
#' 	get_term_descendancy_matrix(hpo, c("HP:0001873", "HP:0001872")),
#' 	matrix(c("HP:0001873","HP:0001872","HP:0001873","HP:0001872"),2,2)
#' )
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
as_row_leaves <- function(term_descendancy_matrix, terms_matrix) {
	terms <- unique(as.vector(terms_matrix))
	term.num.matrix <- apply(
		terms_matrix,
		2,
		function(x) as.integer(factor(x, levels=terms))
	) - 1

	.Call(
		"leaf_matrix",
		term_descendancy_matrix[terms,terms,drop=FALSE],
		term.num.matrix,
		PACKAGE="SimReg"
	)
}

#' @importFrom Rcpp evalCpp
#' @useDynLib SimReg
.log_odds_trace <- function(
	tsm,
	row_is_column_anc,
	as_probs=FALSE,
	x=NULL,
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(row_is_column_anc)))-1,
	g,
	phi,
	gamma,
	alpha_star,
	alpha,
	log_beta,
	logit_mean_f,
	log_alpha_plus_beta_f,
	logit_mean_g,
	log_alpha_plus_beta_g
) {
	result <- .Call(
		"R_log_odds_trace",
		PACKAGE="SimReg",
		tsm,
		row_is_column_anc,
		term_ids,
		case_ids,
		g,
		phi,
		gamma,
		alpha_star,
		alpha,
		log_beta,
		logit_mean_f,
		log_alpha_plus_beta_f,
		logit_mean_g,
		log_alpha_plus_beta_g
	)

	if (as_probs)
		result <- 1 - 1/(1+exp(result))

	result
}	

#' Get trace of log odds of observing the rare genotype (y = 1) for individual cases from the output of \code{\link{sim_reg}}
#'
#' @param term_sim_mat Numeric matrix of similarities between individual terms, typically created with \code{get_term_sim_mat}
#' @template term_descendancy_matrix
#' @template term_sets
#' @param samples Object of class \code{sim_reg_samples}, i.e. the output of the \code{\link{sim_reg}} function.
#' @param g Genotype log odds offsets per individual.
#' @param as_probs Boolean value indicating whether to convert the log odds to probabilities.
#' @return Numeric matrix of log odds trace per individual.
#' @export
log_odds_trace <- function(
	term_sim_mat, 
	term_descendancy_matrix, 
	term_sets, 
	samples, 
	g=rep(FALSE, length(term_sets)), 
	as_probs=FALSE
) {
	stopifnot(class(samples) == "sim_reg_samples")
	.log_odds_trace(
		as_probs=as_probs,
		row_is_column_anc=term_descendancy_matrix,
		x=term_sets,
		tsm=term_sim_mat,
		g=g,
		gamma=samples$gamma,
		alpha=samples$alpha,
		alpha_star=samples$alpha_star,
		log_beta=samples$log_beta,
		logit_mean_f=samples$logit_mean_f,
		log_alpha_plus_beta_f=samples$log_alpha_plus_beta_f,
		logit_mean_g=samples$logit_mean_g,
		log_alpha_plus_beta_g=samples$log_alpha_plus_beta_g,
		phi=t(apply(samples$phi_vector, 1, function(x) match(x, colnames(term_descendancy_matrix))-1))
	)
}

#' Get probabilities y_i given ontological term sets x_i (\code{term_sets}) based on SimReg parameter values sampled with \code{\link{sim_reg}}
#'
#' @param term_sim_mat Numeric matrix of similarities between individual terms, typically created with \code{get_term_sim_mat}
#' @template term_sets
#' @param samples Object of class \code{sim_reg_samples}, i.e. the output of the \code{\link{sim_reg}} function.
#' @param g Genotype log odds offsets per individual.
#' @template ontology
#' @template term_descendancy_matrix
#' @return Numeric vector of probabilities
#' @export
P_y_given_x <- function(
	term_sim_mat,
	term_sets,
	samples,
	g=rep(FALSE, length(term_sets)),
	ontology=NULL,
	term_descendancy_matrix=NULL
) {
	if (is.null(ontology) & is.null(term_descendancy_matrix)) stop("Must supply either 'term_descendancy_matrix' or 'ontology'")
	stopifnot(class(samples) == "sim_reg_samples")
	apply(log_odds_trace(term_sim_mat, term_descendancy_matrix, term_sets, samples, g, TRUE), 2, mean)
}

#' Evaluate asymmetric similarity function
#'
#' @param term_sim_mat Numeric matrix of similarities between individual terms, typically created with \code{get_term_sim_mat}
#' @template term_descendancy_matrix
#' @param phi Character vector of term IDs for phi (or character matrix, where each row is an instance of phi)
#' @template x
#' @template case_ids
#' @template term_ids
#' @param num_cases Number of individuals (automatically determined if \code{x} is given
#' @param average_across_phi Logical value determining whether to use the s_phi function, averaging over similarities to terms in phi, or the s_x, averaging across the best matches to the terms in \code{x}.
#' @return Numeric vector of similarities between phi and \code{x}, or if phi is a matrix, numeric matrix where rows correspond to similarities between x and rows of phi.
#' @export
s <- function(
	term_sim_mat, 
	term_descendancy_matrix, 
	phi,
	x=NULL,
	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length))),
	term_ids=as.integer(match(unlist(x), colnames(term_sim_mat)))-1,
	num_cases=length(x),
	average_across_phi=TRUE
) {
	result <- .Call(
		"R_asym_sim_func",
		PACKAGE="SimReg",
		term_sim_mat,
		term_descendancy_matrix,
		num_cases,
		term_ids,
		case_ids,
		matrix(match(phi, colnames(term_sim_mat))-1, nrow=if (is.matrix(phi)) nrow(phi) else 1, ncol=if (is.matrix(phi)) ncol(phi) else length(phi)),
		average_across_phi	
	)

	if (is.matrix(phi)) result else as.numeric(result)
}

#' Evaluate function s_x
#'
#' @param ... Arguments to pass to \code{\link{s}}
#' @return Similarities
#' @seealso \code{\link{s}}
#' @export
s_x <- function(...) s(average_across_phi=FALSE, ...)

#' Evaluate function s_phi
#'
#' @param ... Arguments to pass to \code{\link{s}}
#' @return Similarities
#' @seealso \code{\link{s}}
#' @export
s_phi <- function(...) s(average_across_phi=TRUE, ...)

#' Extract the estimated posterior probability of gamma = 1
#'
#' @param x A `sim_reg_samples` or `sim_reg_summary` object
#' @return Numeric value
#' @export
`P(gamma=1)` <- function(x) {
	stopifnot(class(x) %in% c("sim_reg_samples", "sim_reg_summary"))
	x[["mean_posterior_gamma"]]
}

#' Gets summarised output of multiple \code{sim_reg} chains
#' 
#' Simple interface to the \code{sim_reg} procedure returning a \code{sim_reg_summary}.
#' 
#' @template sim_reg_description
#' @template ontology
#' @param y Logical vector of genotypes (typically TRUE for rare genotype, FALSE for common genotype).
#' @template x
#' @param g Genotype log odds offset per individual.
#' @param its Number of update cycles to perform .
#' @param chains Number of independent Markov chains used to obtain estimates.
#' @param cores Number of cores to use in parallel computation.
#' @template information_content
#' @template term_descendancy_matrix
#' @template term_sim_mat
#' @param ... Other arguments to be passed to \code{\link{sim_reg}}.
#' @return \code{sim_reg_samples} object.
#' @seealso \code{\link{sim_reg}}
#' @export
#' @importFrom parallel mclapply
sim_reg_mc <- function(
	ontology,
	y,
	x, 
	g=rep(0, length(y)),
	its=10000,
	chains=3,
	cores=1,
	information_content=get_term_info_content(ontology, term_sets=x),
	term_descendancy_matrix=get_term_descendancy_matrix(ontology, names(information_content)),
	term_sim_mat=prune_sim_mat(ontology, get_term_sim_mat(ontology, information_content, term_descendancy_matrix=term_descendancy_matrix)),
	...
) {

	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length)))
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1

	sr <- function(n) sim_reg(
		ontology=ontology,
		y=y,
		x=x,
		g=g,
		its=its,
		information_content=information_content,
		term_descendancy_matrix=term_descendancy_matrix,
		term_sim_mat=term_sim_mat,
		phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
		...
	)
	
	reps <- if (cores > 1 & chains > 1) mclapply(1:chains, mc.cores=cores, FUN=sr) else lapply(1:chains, FUN=sr)

	summary(structure(
		class="sim_reg_samples",
		c(
			mapply(SIMPLIFY=FALSE, FUN=function(tr.name, is.mat) do.call(what=if (is.mat) rbind else c, lapply(reps, "[[", tr.name)), setNames(nm=sim_reg_all_traces), sapply(reps[[1]][sim_reg_all_traces], is.matrix)),
			list(
				mean_posterior_gamma=mean(sapply(reps, "[[", "mean_posterior_gamma")),
				priors=reps[[1]][["priors"]],
				proposal_sds=apply(do.call(what=cbind, lapply(reps, function(rep) simplify2array(rep[["proposal_sds"]]))), 1, mean)
			)
		)
	))
}

