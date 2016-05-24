#' Get logical descendancy matrix for set of terms
#' 
#' @template ontology
#' @template terms
#' @param rows Rows for resultant matrix (defaults to \code{terms})
#' @param cols Cols for resultant matrix (defaults to \code{terms})
#' @return A logical descendancy matrix of \code{terms} by \code{terms} based on DAG structure of ontology. The row term is an ancestor of the column term if result[row.term,col.term] == TRUE
#' @examples 
#' data(hpo)
#' get_term_descendancy_matrix(hpo, c("HP:0001873", "HP:0011877"))
#' @export
#' @importFrom stats setNames
get_term_descendancy_matrix <- function(ontology, terms=NULL, rows=terms, cols=terms) {
	# 'row is column ancestor'
	if (length(terms)==1)
		matrix(FALSE, 1, 1, dimnames=rep(list(terms), 2))
	else
		sapply(
			setNames(cols, cols),
			function(term) setNames(
				rows %in% setdiff(ontology$ancestors[[term]], term),
				rows
			)
		)
}

#' Intersect set of terms with branches of ontology
#' 
#' @template ontology
#' @param branch_roots Character vector of roots of branches you wish to intersect with
#' @template terms
#' @return Character vector of terms
#' @examples 
#' data(hpo)
#' intersection_with_branches(hpo, "HP:0001872", c("HP:0001873", "HP:0011877"))
#' @export
intersection_with_branches <- function(ontology, branch_roots, terms) "if"(
	length(terms) > 0,
	terms[
		sapply(
			ontology$ancestors[terms], 
			function(ancs) any(ancs %in% branch_roots)
		)
	],
	character(0)
)

#' Remove redundant/implied terms
#'
#' @template ontology
#' @template terms
#' @return Character vector of terms
#' @examples
#' data(hpo)
#' minimal_set(hpo, c("HP:0001873", "HP:0001872"))
#' @export
minimal_set <- function(ontology, terms) {
	redundant <- unlist(
		lapply(
			terms,
			function(x) setdiff(ontology$ancestors[[x]], x)
		)
	)
	setdiff(terms,redundant)
}

#' Get set of all ancestors of set of terms
#'
#' @template ontology
#' @template terms
#' @return Character vector of all terms which are an ancestor of at least one term in \code{terms}, including the terms themselves
#' @seealso \code{link{get_descendants}}
#' @examples
#' data(hpo)
#' get_ancestors(hpo, c("HP:0001873", "HP:0011877"))
#' @export
get_ancestors <- function(ontology, terms) {
	unique(
		unlist(
			lapply(terms, function(x) ontology$ancestors[[x]])
		)
	)
}

#' Get set of all descendants of single term
#'
#' @template ontology
#' @param ancestor Character vector of length 1 - the ID of the term whose descendants you wish to retrieve 
#' @param remove_ancestor Boolean indicating whether to remove the given ancestor or not
#' @return Character vector of terms
#' @seealso \code{link{get_ancestors}}
#' @examples
#' data(hpo)
#' get_descendants(hpo, ancestor=c("HP:0001873"))
#' @export
get_descendants <- function(ontology, ancestor, remove_ancestor=FALSE) {
	descs <- unique(
		c(
			ancestor,
			do.call(
				c,
				lapply(ontology$children[[ancestor]], function(child.term) get_descendants(ontology, child.term))
			)
		)
	)

	if (remove_ancestor) setdiff(descs, ancestor) else descs
}

#' Get a matrix with columns of hpo terms and rows of patients,
#'
#' @param term_sets List of character vectors of terms. Result includes only terms which are explicitly present in the list items, so if you wish the result to include even terms which are implicitly present, lapply \code{\link{get_ancestors}} to the argument before passing it to this function
#' @param columns Force result to have these exact columns, entering F for terms which aren't present
#' @return Logical matrix - entry for a patient/hpo term = T if the patient has the term and F otherwise. 
#' @examples
#' get_terms_by_set_matrix(list(Patient1=c("HP:0001873")))
#' @export
get_terms_by_set_matrix <- function(term_sets, columns=NULL) {
	all.terms <- unique(unlist(term_sets))
	result <- t(sapply(term_sets, function(x) all.terms %in% x))
	colnames(result) <- all.terms
	if (!is.null(columns)) {
		result <- cbind(
			result[,intersect(colnames(result), columns)],
			matrix(
				FALSE, 
				nrow(result), 
				length(setdiff(columns, colnames(result))), 
				dimnames=list(rownames(result), setdiff(columns, colnames(result)))
			)
		)[,columns]
	}
	rownames(result) <- names(term_sets)
	result
}

#' Exclude terms descending from particular term from a character vector of terms
#'
#' @template ontology
#' @param branch_root Term whose descendants should be excluded
#' @template terms
#' @return Character vector of terms
#' @export
exclude_branch <- function(ontology, branch_root, terms) Filter(
	x=terms,
	f=function(term) !any(ontology$ancestors[[term]] %in% branch_root)
)

#' Prune all terms descending from given term down to that term and ensure no degeneracy
#'
#' Warning! Prunes down whole branch, even if terms have parents which don't descend from prune point.
#'
#' @template ontology
#' @param prune_to_point Character value equal to term ID which can be included, but whose descendants should be excluded
#' @template terms
#' @return Character vector of terms
#' @export
prune_branch <- function(ontology, prune_to_point, terms) {
	excluded <- exclude_branch(ontology, prune_to_point, terms)
	"if"(
		length(excluded) == length(terms),
		terms,
		c(excluded, prune_to_point)
	)
}

#' Get MPO to HPO list
#'
#' @template ontology
#' @param cross_species_file cross species .obo file, available via http://human-phenotype-ontology.org
#' @return List of HPO terms per MPO term
#' @export
#' @importFrom stats setNames
get_mpo_to_hpo <- function(ontology, cross_species_file) {
	csp.lines <- readLines(cross_species_file)
	term.lines <- which(csp.lines == "[Term]") 

	groupings <- lapply(
		list(
			mpo="^(alt_id|id): (MP:\\d+)",
			hpo="^(alt_id|id|is_a): (HP:\\d+)"
		),
		function(pattern) {
			lines <- grep(pattern, csp.lines)
			split(
				sub(	
					pattern,
					"\\2",
					csp.lines[lines]
				),
				cut(
					lines,
					breaks=c(term.lines, length(csp.lines)+1)
				)
			)
		}
	)

	non.empty <- sapply(
		1:length(groupings$mpo),
		function(x) "&"(
			length(groupings$mpo[[x]]) > 0,
			length(groupings$hpo[[x]]) > 0
		)
	)

	lapply(do.call(what=c, mapply(
		function(mouse.ids, human.ids) setNames(
			rep(list(human.ids), length(mouse.ids)),
			mouse.ids
		),
		groupings$mpo[non.empty],
		groupings$hpo[non.empty],
		USE.NAMES=FALSE
	)), function(x) minimal_set(
			ontology, 
			swap_out_alt_ids(ontology, substr(x, 1, 10))
		)
	)
}

