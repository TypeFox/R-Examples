#' Get human readable, shortened (where possible) HPO term names
#'
#' @template hpo.terms
#' @template terms
#' @return Character vector
#' @examples
#' data(hpo.terms)
#' get.shortened.names(hpo.terms, c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.shortened.names <- function(hpo.terms, terms) gsub(
	"Impaired |(Abnormality of (the )?)|(Abnormal )", 
	"", 
	hpo.terms$name[terms]
) %>% 
(function (x) gsub("^\\s+|\\s+$", "", x)) %>%
sapply(simpleCap)

str.ancs.from.pars <- function(pars, chld) {
	stopifnot(identical(names(pars), names(chld)))
	int.pars <- c(split(as.integer(factor(unlist(pars), levels=names(pars))), unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(pars), sapply(pars, length)))), setNames(nm=setdiff(names(pars), unlist(pars)), rep(list(integer(0)), length(setdiff(names(pars), unlist(pars))))))[names(pars)]
	int.chld <- c(split(as.integer(factor(unlist(chld), levels=names(chld))), unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(chld), sapply(chld, length)))), setNames(nm=setdiff(names(chld), unlist(chld)), rep(list(integer(0)), length(setdiff(names(chld), unlist(chld))))))[names(chld)]

	setNames(nm=names(pars), lapply(ancs.from.pars(
		int.pars,
		int.chld
	), function(x) names(pars)[x]))
}

ancs.from.pars <- function(pars, chld) {
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
#' @import magrittr
get.ontology <- function(file, qualifier="HP") {
	hpo.obo.lines <- readLines(file)

	hpo.term.id.pattern <- paste("^id: (", qualifier, ":\\d+)", sep="")
	hpo.term.name.pattern <- paste("^name: (\\.*)", sep="")
	hpo.term.parent.pattern <- paste("^is_a: (", qualifier, ":\\d+)", sep="")
	hpo.term.pattern <- paste("", qualifier, ":\\d+", sep="")
	hpo.term.alt.id.pattern <- paste("^alt_id: (", qualifier, ":\\d+)", sep="")

	hpo.term.id.lines <- grep(hpo.term.id.pattern, hpo.obo.lines)
	hpo.term.name.lines <- grep(hpo.term.name.pattern, hpo.obo.lines)
	hpo.term.parent.lines <- grep(hpo.term.parent.pattern, hpo.obo.lines)
	hpo.term.alt.id.lines <- grep(hpo.term.alt.id.pattern, hpo.obo.lines)

	hpo.terms <- NULL

	hpo.terms$id <- sub(
		hpo.term.id.pattern,
		"\\1",
		hpo.obo.lines[hpo.term.id.lines]
	)

	hpo.terms$name <- sub(
		hpo.term.name.pattern, 
		"\\1", 
		hpo.obo.lines[hpo.term.name.lines]
	)

	names(hpo.terms$name) <- hpo.terms$id

	Encoding(hpo.terms$name) <- "latin1"
	hpo.terms$name <- iconv(
		hpo.terms$name,
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

	hpo.terms$parents <- split(
		hpo.parent.terms,
		cut(
			hpo.term.parent.lines,
			breaks=c(hpo.term.id.lines,length(hpo.obo.lines)+1),
			labels=hpo.terms$id
		)
	)[hpo.terms$id]
	
	hpo.alt.id.matches <- regexec(
		hpo.term.pattern,
		hpo.obo.lines[hpo.term.alt.id.lines]
	)
	
	hpo.alt.ids <- regmatches(
		hpo.obo.lines[hpo.term.alt.id.lines],
		hpo.alt.id.matches
	)

	hpo.alt.ids <- sapply(hpo.alt.ids, function(x) x[1])
	
	hpo.terms$alt.id <- as.character( 
		cut(
			hpo.term.alt.id.lines,
			breaks=c(hpo.term.id.lines,length(hpo.obo.lines)+1),
			labels=hpo.terms$id
		)
	)	
	names(hpo.terms$alt.id) <- hpo.alt.ids

	names(hpo.terms$id) <- hpo.terms$id

	hpo.terms$children <- c(
		lapply(FUN=as.character, X=split(
			unlist(mapply(SIMPLIFY=FALSE, FUN=rep, names(hpo.terms$parents), sapply(hpo.terms$parents, length))),
			unlist(hpo.terms$parents)
		)),
		setNames(nm=setdiff(hpo.terms$id, unlist(hpo.terms$parents)), rep(list(character(0)), length(setdiff(hpo.terms$id, unlist(hpo.terms$parents)))))
	)[hpo.terms$id]

	hpo.terms$ancestors <- str.ancs.from.pars(hpo.terms$parents, hpo.terms$children)

	hpo.terms$version.info <- head(n=11, hpo.obo.lines)

	hpo.terms
}

#' Remove alternate/deprecated HPO term IDs and swap for new ones
#'
#' @template hpo.terms
#' @template terms
#' @param remove.dead Boolean to indicate whether to strip out terms which can't be found in the given hpo.terms database argument
#' @return A directed adjacency matrix of \code{terms} based on DAG structure of HPO, whereby each term is considered adjacent to it's MRCA in \code{terms}
#' @examples
#' data(hpo.terms)
#' swap.out.alt.ids(hpo.terms, c("HP:0001873"))
#' @export
#' @import magrittr
swap.out.alt.ids <- function(hpo.terms, terms, remove.dead=FALSE) {
	need.swap <- terms %in% names(hpo.terms$alt.id)
	terms[need.swap] <- hpo.terms$alt.id[terms[need.swap]]
	if (remove.dead) 
		terms <- terms[terms %in% hpo.terms$id]
	terms
}

#' Get list of character vector of HPO terms, given character vector of comma separated terms
#'
#' @param character.vector Character vector of comma separated terms
#' @return List of character vectors of HPO terms
#' @examples
#' term.set.list.from.character(c("HP:0001873", "HP:0001902"))
#' @export
#' @import magrittr
term.set.list.from.character <- function(character.vector) strsplit(
	gsub(" ", "", character.vector), 
	split=","
)

#' Get frequency of each term in a set of phenotypes
#'
#' @template hpo.terms
#' @template hpo.phenotypes
#' @param patch.missing Logical indicating whether to include all HPO even if they're not present in the \code{hpo.phenotypes} as if they had occurred once
#' @return Numeric vector of information contents, named by corresponding HPO terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
#' @seealso \code{\link{get.term.info.content}}
#' @examples
#' data(hpo.terms)
#' get.term.frequencies(hpo.terms, list("HP:0001873"))
#' @export
#' @import magrittr
get.term.frequencies <- function(
	hpo.terms, 
	hpo.phenotypes,
	patch.missing=FALSE
) {
	exp(-get.term.info.content(hpo.terms, hpo.phenotypes, patch.missing=FALSE))
}

#' Get information content of each term in a set of phenotypes
#'
#' @template hpo.terms
#' @template hpo.phenotypes
#' @param patch.missing Logical indicating whether to include all HPO even if they're not present in the \code{hpo.phenotypes} as if they had occurred once
#' @return Numeric vector of information contents, named by corresponding HPO terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
#' @examples
#' data(hpo.terms)
#' get.term.info.content(hpo.terms, list("HP:0001873"))
#' @export
#' @import magrittr
get.term.info.content <- function(
	hpo.terms, 
	hpo.phenotypes,
	patch.missing=FALSE
) {
	terms.tab <- table(unlist(lapply(hpo.phenotypes, function(x) get.ancestors(hpo.terms, x))))
	total.patients <- length(hpo.phenotypes)
	terms.numeric <- as.numeric(terms.tab)
	names(terms.numeric) <- names(terms.tab)

	result <- log(total.patients) - ifelse(terms.numeric==0, log(total.patients), log(terms.numeric))
	
	if (patch.missing) {
		#include missing terms and give max information content...
		missing.terms <- setdiff(hpo.terms$id, names(result))
		missing.infos <- rep(max(result), length(missing.terms))
		names(missing.infos) <- missing.terms
		result <- c(result, missing.infos) 
	}

	result
}

#' Object comprising \code{list} of properties of the HPO, indexed by term ID
#' 
#' @name hpo.terms
#' @title HPO Terms object (based on version 887 of the HPO)
#' @docType data
#' @format List of indices containing metadata and structure of HPO 
NULL

#' Object comprising \code{list} of properties of the MPO, indexed by term ID
#' 
#' @name mpo.terms
#' @title MPO Terms object
#' @docType data
#' @format List of indices containing metadata and structure of MPO 
NULL

#' List containing cross-species ontology (MPO to HPO) information - character vectors of HPO terms indexed by associated MPO term IDs
#'
#' @name mpo.to.hpo
#' @title Object containing data for mapping between MPO and HPO
#' @docType data
#' @format List of HPO terms per MPO term
NULL

