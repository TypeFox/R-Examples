#' @importFrom stats setNames
str_ancs_from_pars <- function(pars, chld) {
	stopifnot(identical(names(pars), names(chld)))
	int.pars <- c(split(as.integer(factor(unlist(pars), levels=names(pars))), unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(pars), sapply(pars, length)))), setNames(nm=setdiff(names(pars), unlist(pars)), rep(list(integer(0)), length(setdiff(names(pars), unlist(pars))))))[names(pars)]
	int.chld <- c(split(as.integer(factor(unlist(chld), levels=names(chld))), unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(chld), sapply(chld, length)))), setNames(nm=setdiff(names(chld), unlist(chld)), rep(list(integer(0)), length(setdiff(names(chld), unlist(chld))))))[names(chld)]

	setNames(nm=names(pars), lapply(ancs_from_pars(
		int.pars,
		int.chld
	), function(x) names(pars)[x]))
}

ancs_from_pars <- function(pars, chld) {
	ancs <- as.list(1:length(pars))
	done <- sapply(pars, function(x) length(x) == 0)
	cands <- which(done)
	new.done <- 1:length(cands)
	while (!all(done)) {
		cands <- unique(unlist(chld[cands[new.done]]))
		new.done <- which(sapply(pars[cands], function(x) all(done[x])))
		done[cands[new.done]] <- TRUE
		ancs[cands[new.done]] <- mapply(SIMPLIFY=FALSE, FUN=c, lapply(cands[new.done], function(x) unique(unlist(ancs[pars[[x]]]))), cands[new.done])
	}
	ancs
}

#' Get R-Object representation of ontology from obo file
#'
#' @param file File path of obo file
#' @param qualifier Character vector - "HP" for HPO, "MP" for MPO, etc.
#' @return R-Object (list) representing ontology
#' @export
#' @importFrom stats setNames
get_ontology <- function(file, qualifier="HP") {
	hpo.obo.lines <- readLines(file)

	hpo.term.id.pattern <- paste("^id: (", qualifier, ":\\d+)", sep="")
	hpo.term.name.pattern <- paste("^name: (\\.*)", sep="")
	hpo.term.parent.pattern <- paste("^is_a: (", qualifier, ":\\d+)", sep="")
	hpo.term.pattern <- paste("", qualifier, ":\\d+", sep="")
	hpo.term.alt_id.pattern <- paste("^alt_id: (", qualifier, ":\\d+)", sep="")

	hpo.term.id.lines <- grep(hpo.term.id.pattern, hpo.obo.lines)
	hpo.term.name.lines <- grep(hpo.term.name.pattern, hpo.obo.lines)
	hpo.term.parent.lines <- grep(hpo.term.parent.pattern, hpo.obo.lines)
	hpo.term.alt_id.lines <- grep(hpo.term.alt_id.pattern, hpo.obo.lines)

	ontology <- NULL

	ontology$id <- sub(
		hpo.term.id.pattern,
		"\\1",
		hpo.obo.lines[hpo.term.id.lines]
	)

	ontology$name <- sub(
		hpo.term.name.pattern, 
		"\\1", 
		hpo.obo.lines[hpo.term.name.lines]
	)

	names(ontology$name) <- ontology$id

	Encoding(ontology$name) <- "latin1"
	ontology$name <- iconv(
		ontology$name,
		"latin1",
		"ASCII",
		sub=""
	)

	hpo.parent.term.matches <- regexpr(
		hpo.term.pattern,
		hpo.obo.lines[hpo.term.parent.lines]
	)

	hpo.parent.terms <- substr(
		hpo.obo.lines[hpo.term.parent.lines],
		hpo.parent.term.matches,
		hpo.parent.term.matches + attr(hpo.parent.term.matches, "match.length")-1
	)

	ontology$parents <- split(
		hpo.parent.terms,
		cut(
			hpo.term.parent.lines,
			breaks=c(hpo.term.id.lines,length(hpo.obo.lines)+1),
			labels=ontology$id
		)
	)[ontology$id]
	
	hpo.alt_id.matches <- regexec(
		hpo.term.pattern,
		hpo.obo.lines[hpo.term.alt_id.lines]
	)
	
	hpo.alt_ids <- regmatches(
		hpo.obo.lines[hpo.term.alt_id.lines],
		hpo.alt_id.matches
	)

	hpo.alt_ids <- sapply(hpo.alt_ids, function(x) x[1])
	
	ontology$alt_id <- as.character( 
		cut(
			hpo.term.alt_id.lines,
			breaks=c(hpo.term.id.lines,length(hpo.obo.lines)+1),
			labels=ontology$id
		)
	)	
	names(ontology$alt_id) <- hpo.alt_ids

	names(ontology$id) <- ontology$id

	ontology$children <- c(
		lapply(FUN=as.character, X=split(
			unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(ontology$parents), sapply(ontology$parents, length))),
			unlist(ontology$parents)
		)),
		setNames(nm=setdiff(ontology$id, unlist(ontology$parents)), rep(list(character(0)), length(setdiff(ontology$id, unlist(ontology$parents)))))
	)[ontology$id]

	ontology$ancestors <- str_ancs_from_pars(ontology$parents, ontology$children)

	ontology$version <- hpo.obo.lines[1:11]

	structure(ontology, class="ontology_index")
}

#' Remove alternate/deprecated term IDs and swap for new ones
#'
#' @template ontology
#' @template terms
#' @param remove_dead Boolean to indicate whether to strip out terms which can't be found in the given ontology database argument
#' @return A directed adjacency matrix of \code{terms} based on DAG structure of ontology, whereby each term is considered adjacent to it's MRCA in \code{terms}
#' @examples
#' data(hpo)
#' swap_out_alt_ids(hpo, c("HP:0001873"))
#' @export
swap_out_alt_ids <- function(ontology, terms, remove_dead=FALSE) {
	need.swap <- terms %in% names(ontology$alt_id)
	terms[need.swap] <- ontology$alt_id[terms[need.swap]]
	if (remove_dead) 
		terms <- terms[terms %in% ontology$id]
	terms
}

#' Get frequency of each term in a set of phenotypes
#'
#' @template ontology
#' @template term_sets
#' @param patch_missing Logical indicating whether to include whole ontology even if they're not present in the \code{term_sets} as if they had occurred once
#' @return Numeric vector of information contents, named by corresponding terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
#' @seealso \code{\link{get_term_info_content}}
#' @examples
#' data(hpo)
#' get_term_frequencies(hpo, list("HP:0001873"))
#' @export
get_term_frequencies <- function(
	ontology, 
	term_sets,
	patch_missing=FALSE
) {
	exp(-get_term_info_content(ontology, term_sets, patch_missing=FALSE))
}

#' Get information content of each term in a set of phenotypes
#'
#' @template ontology
#' @template term_sets
#' @param patch_missing Logical indicating whether to include all ontology terms even if they're not present in the \code{term_sets} as if they had occurred once
#' @return Numeric vector of information contents, named by corresponding terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
#' @examples
#' data(hpo)
#' get_term_info_content(hpo, list("HP:0001873"))
#' @export
get_term_info_content <- function(
	ontology, 
	term_sets,
	patch_missing=FALSE
) {
	terms.tab <- table(unlist(lapply(term_sets, function(x) get_ancestors(ontology, x))))
	total.patients <- length(term_sets)
	terms.numeric <- as.numeric(terms.tab)
	names(terms.numeric) <- names(terms.tab)

	result <- log(total.patients) - ifelse(terms.numeric==0, log(total.patients), log(terms.numeric))
	
	if (patch_missing) {
		#include missing terms and give max information content...
		missing.terms <- setdiff(ontology$id, names(result))
		missing.infos <- rep(max(result), length(missing.terms))
		names(missing.infos) <- missing.terms
		result <- c(result, missing.infos) 
	}

	result
}

#' \code{ontology_index} object encapsulating structure of the Human Phenotype Ontology (HPO) comprising a \code{list} of lists/vectors of properties of HPO terms indexed by term ID
#' 
#' @name hpo
#' @title HPO Terms object
#' @docType data
#' @format List of lists and vectors
NULL

#' \code{ontology_index} object encapsulating structure of the Human Phenotype Ontology (MPO) comprising a \code{list} of lists/vectors of properties of MPO terms indexed by term ID
#' 
#' @name mpo
#' @title MPO Terms object (based on version 887 of the MPO)
#' @docType data
#' @format List of lists and vectors
NULL

#' List containing cross-species ontology (MPO to HPO) information - character vectors of HPO terms indexed by associated MPO term IDs
#'
#' @name mpo_to_hpo
#' @title Object containing data for mapping between MPO and HPO
#' @docType data
#' @format List of HPO terms per MPO term
NULL

