#' Get matrix of similarity rank from similarity matrix
#'
#' Given a lower triangular similarity matrix, construct a distance matrix where the rows are the ranks of the column cases, ranked by similarity to the row case
#'
#' @param similarity_matrix Lower triangular numeric matrix of similarities, where the rownames and colnames are identical to the case IDs
#' @return Rank matrix (not necessarily symmetrical
#' @export
get_similarity_rank_matrix <- function(similarity_matrix) { 
	similarity_matrix[upper.tri(similarity_matrix, diag=FALSE)] <- 0
	similarity_matrix <- similarity_matrix + t(similarity_matrix)
	diag(similarity_matrix) <- NA
	ranked.sim.matrix <- apply(similarity_matrix, 1, function(x) rank(-x, ties.method="max", na.last="keep"))
	t(ranked.sim.matrix)
}

#' Get term-term similarity matrix 
#'
#' @template ontology
#' @template information_content 
#' @template term_descendancy_matrix
#' @param method Character vector - either "lin" (to use Lin's expression of similarity of terms) or "resnik" (to use Resnik's definition)
#' @param restrict_rows_to Restrict rows to Character vector of term names for which to compute corresponding similarities to all terms \code{information_content}
#' @return Numeric matrix of term-term simliarities
#' @seealso \code{\link{resnik}}, \code{\link{lin}}
#' @importFrom ontologyIndex get_ancestors get_term_descendancy_matrix
#' @export
get_term_sim_mat <- function(
	ontology, 
	information_content, 
	term_descendancy_matrix=get_term_descendancy_matrix(ontology, rows=get_ancestors(ontology, restrict_rows_to), cols=names(information_content)), 
	method=c("lin", "resnik"), 
	restrict_rows_to=names(information_content)
) {
	result <- matrix(0, nrow=length(restrict_rows_to), ncol=length(information_content), dimnames=list(restrict_rows_to, names(information_content)))

	for (term in get_ancestors(ontology, restrict_rows_to)) {
		#descs <- c(match(term, names(information_content)), which(term_descendancy_matrix[term,]))
		descs <- c(term, names(which(term_descendancy_matrix[term,])))
		in.rows <- intersect(restrict_rows_to, descs)
		result[in.rows,descs] <- pmax(matrix(information_content[term], nrow=length(in.rows), ncol=length(descs)), result[in.rows,descs])
	}

	if (method[1] == "lin") {
		ic.sum <- outer(information_content[restrict_rows_to], information_content, "+")

		result <- 2 * result / ic.sum

		result[which(arr.ind=TRUE, is.na(result))] <- 0
	}

	result
}

#' Get similarity p-value
#'
#' @template sym_sim_mat
#' @template case_group
#' @template min_its 
#' @template max_its
#' @template signif
#' @template log_dismiss
#' @return Numeric p-value of similarity computed from null distribution (group similarity of all groups of same size as \code{case_group}
#' @seealso \code{\link{get_stratified_sim_p}}
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_p <- function(sym_sim_mat, case_group, min_its=1e3, max_its=1e6, signif=0.05, log_dismiss=-15) {
	.Call(
		"sim_p",
		sym_sim_mat,
		"if"(class(case_group) == "character", which(rownames(sym_sim_mat) %in% case_group)-1, case_group-1),
		min_its,
		max_its,
		signif,
		log_dismiss,
		PACKAGE="ontologySimilarity"
	)
}

#' Get stratified similarity p-value
#'
#' @template sym_sim_mat
#' @param row_strata_factor Factor of strata to which the rows of the similarity matrix belong
#' @template case_group
#' @param min_its Minimum number of simulated group similarities to calculate
#' @param max_its Maximum number of simulated group similarities to calculate
#' @template signif
#' @template log_dismiss
#' @param subgroup_size Calculate significance of similarity of most similar subgroup within \code{case_group}. Set to NULL to use the whole group (default).
#' @return Numeric p-value of similarity computed from null distribution (group similarity of all groups containing the same number of cases from each strata)
#' @export
get_stratified_sim_p <- function(
	sym_sim_mat,
	row_strata_factor=factor(rep(1, nrow(sym_sim_mat))),
	case_group,
	min_its=1e3,
	max_its=1e6,
	signif=0.05,
	log_dismiss=-15,
	subgroup_size=NULL
 ) {
	stopifnot(class(row_strata_factor) == "factor")
	char.group <- "if"(class(case_group) != "character", rownames(sym_sim_mat)[case_group], case_group)
	new.ordering <- order(row_strata_factor)
	get_stratified_sim_p_sorted(
		sym_sim_mat[new.ordering, new.ordering],
		row_strata_factor[new.ordering],
		which(rownames(sym_sim_mat)[new.ordering] %in% char.group),
		min_its,
		max_its,
		signif,
		log_dismiss,
		subgroup_size
	)
}

#' Get stratified similarity p-value (sorted)
#'
#' @param sorted_sym_sim_mat Symmetric numeric term set similarity matrix with rows sorted by strata
#' @param row_strata_factor Factor of strata to which the rows of the similarity matrix belong
#' @template case_group
#' @param min_its Minimum number of simulated group similarities to calculate
#' @param max_its Maximum number of simulated group similarities to calculate
#' @template signif
#' @template log_dismiss
#' @param subgroup_size Calculate significance of similarity of most similar subgroup within \code{case_group}. Set to NULL to use the whole group (default).
#' @return Numeric p-value of similarity computed from null distribution (group similarity of all groups containing the same number of cases from each strata)
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_stratified_sim_p_sorted <- function(
	sorted_sym_sim_mat,
	row_strata_factor,
	case_group,
	min_its=1e3,
	max_its=1e6,
	signif=0.05,
	log_dismiss=-15,
	subgroup_size=NULL
) {
	stopifnot(class(case_group) == "integer")

	"if"(
		is.null(subgroup_size),
		.Call(
			"stratified_sim_p",
			sorted_sym_sim_mat,
			case_group-1,
			as.integer(table(row_strata_factor)),
			as.integer(table(row_strata_factor[case_group])),
			min_its,
			max_its,
			signif,
			log_dismiss,
			PACKAGE="ontologySimilarity"
		),
		.Call(
			"stratified_best_subgroup_p",
			subgroup_size,
			sorted_sym_sim_mat,
			case_group-1,
			as.integer(table(row_strata_factor)),
			as.integer(table(row_strata_factor[case_group])),
			min_its,
			max_its,
			signif,
			log_dismiss,
			PACKAGE="ontologySimilarity"
		)
	)
}

#' Get `term sets to term' similarity matrix
#'
#' Using a matrix of between-term similarities (e.g. the kind obtained from applying the function \code{\link{get_term_sim_mat}}), create a numeric matrix of `term set-to-term' similarities using either the `best-match' approach (i.e. the similarity of the term will be the maximum term-to-term similarity in \code{term_sim_mat} amongst the term set terms).
#'
#' @template term_sim_mat
#' @template term_sets
#' @return Numeric matrix of term set-to-term similarities
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_term_set_to_term_sims <- function(
	term_sim_mat,
	term_sets
) {

	term.ids=as.integer(match(unlist(term_sets), colnames(term_sim_mat)))-1
	case.ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(term_sets)-1), sapply(term_sets, length)))
	num_cases=length(term_sets)
	phenotype.labels=names(term_sets)

	result <- t(.Call(
		"R_get_sim_matrix_2_sets",
		PACKAGE="ontologySimilarity",
		(1:ncol(term_sim_mat))-1,
		(1:ncol(term_sim_mat))-1,
		ncol(term_sim_mat),
		term.ids,
		case.ids,
		num_cases,
		term_sim_mat
	))

	rownames(result) <- phenotype.labels
	colnames(result) <- colnames(term_sim_mat)

	result
}

#' Get similarity matrix from a term-to-term similarity matrix between two sets of term sets
#'
#' Equivalent to \code{\link{get_sim_mat}} for two sets of term sets. 
#'
#' @template term_sim_mat
#' @template term_sets
#' @param term_sets2 Second set of term sets
#' @param combine Character string - either \code{average} or \code{product}, indicating whether to use the `best-match-average' or `best-match-product' method
#' @return Numeric matrix of between-term set similarities
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_grid <- function(
	term_sim_mat,
	term_sets,
	term_sets2=term_sets,
	combine=c("average", "product")
) {
	term.ids=as.integer(match(unlist(term_sets), colnames(term_sim_mat)))-1
	case.ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(term_sets)-1), sapply(term_sets, length)))
	num.cases=length(term_sets)
	phenotype.labels=names(term_sets)

	term.ids2=as.integer(match(unlist(term_sets2), colnames(term_sim_mat)))-1
	case.ids2=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(term_sets2)-1), sapply(term_sets2, length)))
	num.cases2=length(term_sets2)
	phenotype.labels2=names(term_sets2)

	result <- (if (combine[1] == "average") function(x, y) { (x + y) / 2 } else match.fun("*"))(
		.Call(
			"R_get_sim_matrix_2_sets",
			PACKAGE="ontologySimilarity",
			term.ids,
			case.ids,
			num.cases,
			term.ids2,
			case.ids2,
			num.cases2,
			term_sim_mat
		),
		t(.Call(
			"R_get_sim_matrix_2_sets",
			PACKAGE="ontologySimilarity",
			term.ids2,
			case.ids2,
			num.cases2,
			term.ids,
			case.ids,
			num.cases,
			term_sim_mat
		))
	)

	rownames(result) <- phenotype.labels
	colnames(result) <- phenotype.labels2
	result
}

#' Get similarity matrix from a term similarity matrix
#'
#' Using a matrix of between-term similarities (e.g. the kind obtained from applying the function \code{\link{get_term_sim_mat}}), create a numeric matrix of `between-term set' similarities, using either the `best-match-average' or `best-match-product' approach (i.e. where the 2 scores obtained by applying the asymmetric `best-match' similarity function to two term sets in each order are combined by taking the average or the product respectively).
#'
#' @template term_sim_mat
#' @template term_sets
#' @param combine Character string - either \code{average} or \code{product}, indicating whether to use the `best-match-average' or `best-match-product' method
#' @return Numeric matrix of between-term set similarities
#' @examples
#' suppressPackageStartupMessages(library(ontologyIndex))
#' data(hpo)
#' set.seed(1)
#' #random set of terms with ancestors
#' terms <- get_ancestors(hpo, sample(hpo$id, size=30))
#' #set information content of terms (as if each term occurs with frequency `1/n`)
#' information_content <- get_term_info_content(hpo, term_sets=as.list(terms))
#' #similarity of term pairs
#' tsm <- get_term_sim_mat(hpo, information_content)
#' #5 random term sets (call them *phenotypes*) with (at most) 8 terms (removing redundant ones)
#' phenotypes <- lapply(replicate(simplify=FALSE, n=5, 
#'   expr=sample(terms, size=8)), minimal_set, ontology=hpo)
#' get_sim_mat(tsm, phenotypes)
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_mat <- function(
	term_sim_mat,
	term_sets,
	combine=c("average", "product")
) {

	term.ids=as.integer(match(unlist(term_sets), colnames(term_sim_mat)))-1
	case.ids=unlist(mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(term_sets)-1), sapply(term_sets, length)))
	num.cases=length(term_sets)
	phenotype.labels=names(term_sets)

	un.sym <- .Call(
		"R_get_sim_matrix",
		PACKAGE="ontologySimilarity",
		term.ids,
		case.ids,
		num.cases,
		term_sim_mat
	)

	result <- if (combine[1] == "average") (un.sym + t(un.sym))/2 else if (combine[1] == "product") un.sym * t(un.sym) else stop("Combine method must be average or product!")
	rownames(result) <- colnames(result) <- phenotype.labels
	result
}

#' Calculate the group similarity of a set of row/column indices
#'
#' @template sym_sim_mat
#' @param group Integer/character vector giving the indices/names of rows/columns in the \code{sym_sim_mat} for which to compute a similarity 
#' @return Numeric value of group similarity
#' @export
get_sim <- function(
	sym_sim_mat,
	group
) {
	mean(sym_sim_mat[group,group], na.rm=TRUE)
}

#' Get random sample of group similarities for groups of given size
#'
#' @template sym_sim_mat
#' @param group_size Integer giving the number of members of a group
#' @param sample_size Number of samples to draw
#' @return Numeric vector of random group similarities
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_sample <- function(
	sym_sim_mat,
	group_size,
	sample_size
) {
	.Call(
		"R_get_sim_sample",
		PACKAGE="ontologySimilarity",
		sym_sim_mat,
		group_size,
		sample_size
	)
}

#' Get p-value for best-match-average similarity between two term sets 
#'
#' Get p-value by comparison of similarity of \code{profile} and \code{term_set} to distribution of similarities obtained by simulating random term sets (i.e. comprising random sets of terms the same size as \code{term_set}) and calculating the similarity to \code{profile}.
#'
#' @template term_sim_mat
#' @param profile Character vector of ontological term IDs
#' @param term_set Character vector of ontological term IDs
#' @param sample_from Character vector of ontological term IDs to sample term sets from
#' @template min_its 
#' @template max_its
#' @template signif
#' @template log_dismiss
#' @return Numeric p-value
#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_to_profile_p <- function(
	term_sim_mat,
	profile,
	term_set,
	sample_from=colnames(term_sim_mat),
	min_its=1e3,
	max_its=1e5,
	signif=0.05,
	log_dismiss=-10
) {
	.Call(
		'R_get_sim_to_profile_p',
		term_sim_mat,
		match(sample_from, colnames(term_sim_mat))-1,
		match(profile, rownames(term_sim_mat))-1,
		match(term_set, colnames(term_sim_mat))-1,
		min_its,
		max_its,
		signif,
		log_dismiss
	)
}

resnik_asym <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	mean(sapply(
		term_set_1,	
		function(x) max(information_content[
			intersect(ontology$ancestors[[x]], get_ancestors(ontology, term_set_2))
		])
	))
}

#' Calculate Resnik similarity score of two term sets
#'
#' Warning! This function is slow - performing large numbers of `between term-set' similarity calculations should be done using \code{\link{get_sim_grid}}.
#'
#' @template ontology
#' @template information_content
#' @param term_set_1 Character vector of terms
#' @param term_set_2 Character vector of terms
#' @seealso \code{\link{lin}}, \code{\link{get_term_sim_mat}}
#' @references Resnik, P. (1995). Using information content to evaluate semantic similarity in a taxonomy. Proceedings of the 14th IJCAI 1, 448-453.
#' @export
resnik <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	(resnik_asym(ontology, information_content, term_set_1, term_set_2) + resnik_asym(ontology, information_content, term_set_2, term_set_1))/2	
}

lin_asym <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	mean(sapply(
		term_set_1,	
		function(t1) {
			max(sapply(term_set_2, function(t2) 2*max(information_content[intersect(ontology$ancestors[[t2]], ontology$ancestors[[t1]])])/(information_content[t1]+information_content[t2])))
		}
	))
}

#' Calculate Lin similarity score of two term sets
#'
#' Warning! This function is slow - performing large numbers of `between term-set' similarity calculations should be done using \code{\link{get_sim_grid}}.
#'
#' @template ontology
#' @template information_content
#' @param term_set_1 Character vector of terms
#' @param term_set_2 Character vector of terms
#' @seealso \code{\link{resnik}}, \code{\link{get_term_sim_mat}}
#' @references Lin, D. (1998). An information-theoretic definition of similarity. In Shavlik, J. W., ed., Proceedings of the Fifteenth International Conference on Machine Learning (ICML 1998), Madison, Wisconsin, USA, July 24-27, 1998. (Morgan Kaufmann) pp. 296-304. 
#' @export
lin <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	(lin_asym(ontology, information_content, term_set_1, term_set_2) + lin_asym(ontology, information_content, term_set_2, term_set_1))/2	
}
