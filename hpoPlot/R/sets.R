#' Get an adjacency to MRCA matrix for set of HPO terms
#' 
#' @template hpo.terms
#' @template terms
#' @return A logical matrix represenging the directed adjacency matrix of \code{terms} based on DAG structure of HPO, whereby a TRUE entry signifies the term correspnding to the column is MRCA of the row term in \code{terms}
#' @seealso \code{\link{get.term.adjacency.matrix}}
#' @examples
#' data(hpo.terms)
#' get.term.pseudo.adjacency.matrix(hpo.terms, c("HP:0000118", "HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.term.pseudo.adjacency.matrix <- function(hpo.terms, terms) setDimNames(
	t(
		sapply(
			setNames(terms, terms),
			function(term) "%in%"(
				setNames(terms, terms),
				clean.terms(	
					hpo.terms, 
					setdiff(
						intersect(terms, hpo.terms$ancestors[[term]]),
						term
					)
				)
			)
		)
	),
	rep(
		list(terms),
		2
	)
)

#' Get an adjacency for set of HPO terms
#' 
#' @template hpo.terms
#' @template terms
#' @return A logical matrix representing the directed adjacency matrix of \code{terms} based on DAG structure of HPO, whereby a TRUE entry signifies that the term corresponding to the column is a parent term of the term correspnding to the row.
#' @seealso \code{\link{get.term.pseudo.adjacency.matrix}}
#' @examples
#' data(hpo.terms)
#' get.term.adjacency.matrix(hpo.terms, c("HP:0000118", "HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.term.adjacency.matrix <- function(hpo.terms, terms) {
	names(terms) <- terms
	adj.mat <- sapply(
		terms,
		function(term) terms %in% hpo.terms$parents[[term]]
	)
	rownames(adj.mat) <- colnames(adj.mat)
	adj.mat <- t(adj.mat)
	return(adj.mat)
}

#' Get logical descendancy matrix for set of terms
#' 
#' @template hpo.terms
#' @template terms
#' @param rows Rows for resultant matrix (defaults to \code{terms})
#' @param cols Cols for resultant matrix (defaults to \code{terms})
#' @return A logical descendancy matrix of \code{terms} by \code{terms} based on DAG structure of HPO, where by the row  term is an ancestor of the column term if result[row.term,col.term] == TRUE
#' @examples 
#' data(hpo.terms)
#' get.term.descendancy.matrix(hpo.terms, c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.term.descendancy.matrix <- function(hpo.terms, terms=NULL, rows=terms, cols=terms) {
	# 'row is column ancestor'
	if (length(terms)==1)
		matrix(FALSE, 1, 1, dimnames=rep(list(terms), 2))
	else
		sapply(
			setNames(cols, cols),
			function(term) setNames(
				rows %in% setdiff(hpo.terms$ancestors[[term]], term),
				rows
			)
		)
}

#' Intersect set of terms with branches of HPO
#' 
#' @template hpo.terms
#' @param branch.roots Character vector of roots of branches you wish to intersect with
#' @template terms
#' @return Character vector of terms
#' @examples 
#' data(hpo.terms)
#' intersection.with.branches(hpo.terms, "HP:0001872", c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
intersection.with.branches <- function(hpo.terms, branch.roots, terms) "if"(
	length(terms) > 0,
	terms[
		sapply(
			hpo.terms$ancestors[terms], 
			function(ancs) any(ancs %in% branch.roots)
		)
	],
	character(0)
)

#' Remove redundant/implied terms
#'
#' @template hpo.terms
#' @template terms
#' @return Character vector of HPO terms
#' @examples
#' data(hpo.terms)
#' clean.terms(hpo.terms, c("HP:0001873", "HP:0001872"))
#' @export
#' @import magrittr
clean.terms <- function(hpo.terms, terms) {
	redundant <- unlist(
		lapply(
			terms,
			function(x) setdiff(hpo.terms$ancestors[[x]], x)
		)
	)
	setdiff(terms,redundant)
}

#' Get set of all ancestors of set of terms
#'
#' @template hpo.terms
#' @template terms
#' @return Character vector of all HPO terms which are an ancestor of at least one term in \code{terms}, including the terms themselves
#' @seealso \code{link{get.descendants}}
#' @examples
#' data(hpo.terms)
#' get.ancestors(hpo.terms, c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.ancestors <- function(hpo.terms, terms) {
	unique(
		unlist(
			lapply(terms, function(x) hpo.terms$ancestors[[x]])
		)
	)
}

#' Get set of all descendants of single term
#'
#' @template hpo.terms
#' @param ancestor Character vector of length 1 - the HPO code of the term whose descendants you wish to retrieve 
#' @param remove.ancestor Boolean indicating whether to remove the given ancestor or not
#' @return Character vector of terms
#' @seealso \code{link{get.ancestors}}
#' @examples
#' data(hpo.terms)
#' get.descendants(hpo.terms, ancestor=c("HP:0001873"))
#' @export
#' @import magrittr
get.descendants <- function(hpo.terms, ancestor, remove.ancestor=FALSE) unique(
	c(
		ancestor,
		do.call(
			c,
			lapply(hpo.terms$children[[ancestor]], function(child.term) get.descendants(hpo.terms, child.term))
		)
	)
) %>% 
(function(descs) "if"(
	remove.ancestor,
	setdiff(descs, ancestor),
	descs
))
	
#' Get a matrix with columns of hpo terms and rows of patients,
#'
#' @param hpo.phenotypes List of character vectors of HPO terms. Result includes only terms which are explicitly present in the list items, so if you wish the result to include even terms which are implicitly present, lapply \code{\link{get.ancestors}} to the argument before passing it to this function
#' @param columns Force result to have these exact columns, entering F for terms which aren't present
#' @return Logical matrix - entry for a patient/hpo term = T if the patient has the term and F otherwise. 
#' @examples
#' get.case.term.matrix(list(Patient1=c("HP:0001873")))
#' @export
#' @import magrittr
get.case.term.matrix <- function(hpo.phenotypes, columns=NULL) {
	all.terms <- unique(unlist(hpo.phenotypes))
	result <- t(sapply(hpo.phenotypes, function(x) all.terms %in% x))
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
	result
}

#' Get a minimal set of terms which can be used to partition a set of phenotypes 
#'
#' @template hpo.terms
#' @template hpo.phenotypes
#' @return Character vector of set of terms, excluding terms for which the presence of their descendants all partition the set of terms in the same way
#' @examples
#' data(hpo.terms)
#' remove.uninformative.terms(hpo.terms, list(Patient1=c("HP:0001873")))
#' @export
#' @import magrittr
remove.uninformative.terms <- function(hpo.terms, hpo.phenotypes) {
	with.ancs <- lapply(hpo.phenotypes, function(x) get.ancestors(hpo.terms, x))
	all.terms <- unique(unlist(with.ancs))
	terms <- Filter(
		x=unique(unlist(with.ancs)),
		f=function(term) {
			if (
				"=="(
					length(
						intersect(
							hpo.terms$children[[term]],
							all.terms
						)
					),
					0
				)
			) return(TRUE)

			patients.of.each.child <- lapply(
				intersect(hpo.terms$children[[term]], all.terms),
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
	
	return(terms)
}

#' Exclude terms descending from particular term from a character vector of terms
#'
#' @template hpo.terms
#' @param branch.root HPO term whose descendants should be excluded
#' @template terms
#' @return Character vector of terms
#' @export
#' @import magrittr
exclude.branch <- function(hpo.terms, branch.root, terms) Filter(
	x=terms,
	f=function(term) !any(hpo.terms$ancestors[[term]] %in% branch.root)
)

#warning! prunes down whole branch, even if terms have parents which don't descend from prune point... 
#' Prune all terms descending from given term down to that term and ensure no degeneracy
#'
#' @template hpo.terms
#' @param prune.to.point HPO term which can be included, but whose descendants should be excluded
#' @template terms
#' @return Character vector of terms
#' @export
#' @import magrittr
prune.branch <- function(hpo.terms, prune.to.point, terms) {
	excluded <- exclude.branch(hpo.terms, prune.to.point, terms)
	"if"(
		length(excluded) == length(terms),
		terms,
		c(excluded, prune.to.point)
	)
}

#' Get MPO to HPO R-Object
#'
#' @template hpo.terms
#' @param cross.species.file cross species obo file, downloadable from http://compbio.charite.de/hudson/? website
#' @return List of HPO terms per MPO term
#' @export
#' @import magrittr
get.mpo.to.hpo <- function(hpo.terms, cross.species.file) {
	csp.lines <- readLines(cross.species.file)
	term.lines <- which(csp.lines == "[Term]") 

	groupings <- lapply(
		#patterns
		list(
			mpo="^(alt_id|id): (MP:\\d+)",
			hpo="^(alt_id|id|is_a): (HP:\\d+)"
		),
		function(pattern) 
			grep(pattern, csp.lines) %>%
			(function(lines) split(
				sub(	
					pattern,
					"\\2",
					csp.lines[lines]
				),
				cut(
					lines,
					breaks=c(term.lines, length(csp.lines)+1)
				)
			))
	)

	#ensure that no overlap of mammalian phenotype ids...
	stopifnot(
		"=="(
			Filter(x=groupings$mpo, f=function(x) length(x) > 0) %>% length,
			Filter(x=groupings$mpo, f=function(x) length(x) > 0) %>% 
			c %>% unique %>% length
		)
	)

	non.empty <- sapply(
		1:length(groupings$mpo),
		function(x) "&"(
			length(groupings$mpo[[x]]) > 0,
			length(groupings$hpo[[x]]) > 0
		)
	)

	mapply(
		function(mouse.ids, human.ids) setNames(
			rep(list(human.ids), length(mouse.ids)),
			mouse.ids
		),
		groupings$mpo[non.empty],
		groupings$hpo[non.empty],
		USE.NAMES=FALSE
	) %>% 
	#collapse into single list... safe in the knowledge that all the mouse phenotypes are unique...
	do.call(what=c) %>% 
	#hack to remove trailing parent term descriptions...
	lapply(function(x) substr(x, 1, 10)) %>%
	lapply(
		function(x) clean.terms(
			hpo.terms, 
			swap.out.alt.ids(hpo.terms, x)
		)
	)
}

