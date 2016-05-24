#' Get an adjacency matrix for a set of ontological terms
#' 
#' @template ontology
#' @template terms
#' @return A logical matrix representing the adjacency matrix of \code{terms} based on the directed acyclic graph of \code{ontology}. A \code{TRUE} entry means the term correspnding to the column is a parent of the row term within \code{terms}.
#' @seealso \code{\link{get_adjacency_matrix}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_pseudo_adjacency_matrix(hpo, c("HP:0000118", "HP:0001873", "HP:0011877"))
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
get_pseudo_adjacency_matrix <- function(ontology, terms) structure(
	t(
		sapply(
			setNames(terms, terms),
			function(term) "%in%"(
				setNames(terms, terms),
				minimal_set(	
					ontology, 
					setdiff(
						intersect(terms, ontology$ancestors[[term]]),
						term
					)
				)
			)
		)
	),
	dimnames=rep(
		list(terms),
		2
	)
)

#' Get an adjacency matrix for a set of ontological terms
#' 
#' @template ontology
#' @template terms
#' @return A logical matrix representing the adjacency matrix of \code{terms} based on the directed acyclic graph of \code{ontology}. A \code{TRUE} entry means the term correspnding to the column is a parent of the row term in \code{ontology}.
#' @seealso \code{\link{get_pseudo_adjacency_matrix}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_adjacency_matrix(hpo, c("HP:0000118", "HP:0001873", "HP:0011877"))
#' @export
#' @import ontologyIndex
get_adjacency_matrix <- function(ontology, terms) {
	names(terms) <- terms
	adj.mat <- sapply(
		terms,
		function(term) terms %in% ontology$parents[[term]]
	)
	rownames(adj.mat) <- colnames(adj.mat)
	t(adj.mat)
}

#' Retain only the most specific terms which are present in each unique set of term sets
#'
#' Useful in finding a relatively small set of terms which captures the structure and overlap of terms within a set of term sets.
#'
#' @template ontology
#' @template term_sets
#' @return Character vector of terms
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' remove_uninformative_terms(hpo, list(Patient1=c("HP:0001873")))
#' @export
#' @import ontologyIndex
remove_uninformative_terms <- function(ontology, term_sets) {
	with.ancs <- lapply(term_sets, function(x) get_ancestors(ontology, x))
	all.terms <- unique(unlist(with.ancs))
	terms <- Filter(
		x=unique(unlist(with.ancs)),
		f=function(term) {
			if (
				"=="(
					length(
						intersect(
							ontology$children[[term]],
							all.terms
						)
					),
					0
				)
			) return(TRUE)

			patients.of.each.child <- lapply(
				intersect(ontology$children[[term]], all.terms),
				function(child) sapply(
					with.ancs, 
					function(patient.terms) child %in% patient.terms
				)
			)
			patients.of.term <- sapply(
				with.ancs,
				function(patient.terms) term %in% patient.terms
			)
			!all(
				sapply(
					patients.of.each.child,
					function(x) identical(
						x,
						patients.of.term
					)
				)
			)
		}
	)
	
	terms
}

