#' Similarity regression on multiple \code{y} vectors
#'
#' Applies \code{\link{sim_reg}} regression to multiple binary 'genotypes', \code{y}, encoded in a logical matrix having the individual binary genotypes as rows, against the same set of ontological 'phenotypes'. Thus, the individual subject's genotypes occupy the columns of \code{y} and the number of columns is equal to the length of the list of phenotypes \code{x}. The user can supply a function to be applied to the \code{sim_reg_samples} object for each completed inference. Suitable choices for the function include \code{\link{summary}} (default) and \code{identity} (though identity will lead to extremely large objects if there are many rows in the \code{y} matrix); the number of cores to use.
#'
#' @template ontology
#' @param y Logical matrix of genotypes (typically TRUE for rare genotype, FALSE for common genotype).
#' @template x
#' @param g Genotype log odds offset per individual.
#' @param lit_sim_mat Numeric matrix of prior weights for terms inclusion in phi for each row of \code{y}, where the columns match to those in \code{term_sim_mat}. Thus, must have exactly the same number of rows as \code{y}, and same number of columns as \code{term_sim_mat}. Defaults to 'all equal' (=1).
#' @param summary_function Function to apply to the samples generated conditioning on each row of \code{y}.
#' @param simplify Logical value determining whether result is simplified by using \code{rbind} to turn it into a \code{data.frame}.
#' @param mc_cores Number of cores to use - passed to \code{mcmapply} function.
#' @template information_content
#' @template term_descendancy_matrix
#' @template term_sim_mat
#' @param ... Other arguments to be passed to \code{\link{sim_reg}}.
#' @return List or data frame of summarised output of \code{\link{sim_reg}}.
#' @export
#' @importFrom parallel mcmapply 
#' @importFrom stats rnorm runif
#' @importFrom ontologyIndex get_term_descendancy_matrix
multi_sim_reg <- function(
	ontology,
	x,
	y,
	g=matrix(0, nrow=nrow(y), ncol=ncol(y)),
	lit_sim_mat=NULL,
	summary_function=summary,
	simplify=FALSE,
	mc_cores=1L,
	information_content=get_term_info_content(ontology, term_sets=x),
	term_descendancy_matrix=get_term_descendancy_matrix(ontology, names(information_content)),
	term_sim_mat=prune_sim_mat(ontology, get_term_sim_mat(ontology, information_content, term_descendancy_matrix=term_descendancy_matrix)),
	...
) {
	stopifnot(mode(y) == "logical" & class(y) == "matrix")
	stopifnot(length(x) == ncol(y))
	stopifnot(nrow(y) > 0)

	case_ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(x)-1), sapply(x, length)))
	term_ids=as.integer(match(unlist(x), colnames(term_descendancy_matrix)))-1

	result <- mcmapply(
		mc.cores=mc_cores,
		SIMPLIFY=FALSE,

		FUN=function(y_row, g_row, lit_sims_row) summary_function(
			sim_reg(
				y=y_row,
				g=g_row,
				x=x,
				lit_sims=lit_sims_row,
				information_content=information_content,
				term_descendancy_matrix=term_descendancy_matrix,
				term_sim_mat=term_sim_mat,
				case_ids=case_ids,
				term_ids=term_ids,
				phi_jumps=c(0:(ncol(term_descendancy_matrix)-1), rep(match(unlist(lapply(x[y], get_ancestors, ontology=ontology)), colnames(term_descendancy_matrix))-1, 50)),
				...
			)
		),
		split(y, seq(nrow(y))),
		split(g, seq(nrow(g))),
		if (is.null(lit_sim_mat)) rep(list(NULL), nrow(y)) else split(lit_sim_mat, seq(nrow(lit_sim_mat)))
	)

	if (simplify) do.call(what=rbind, lapply(result, as.data.frame)) else result
}

