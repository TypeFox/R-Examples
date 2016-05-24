#' Select \code{n} most prevalent terms in \code{term_sets}
#'
#' @template ontology
#' @template term_sets
#' @param n Integer 
#' @template terms
#' @return Character vector of length at most \code{n}
#' @seealso \code{\link{remove_terms_with_less_than_n_occurrences}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' n_most_frequent_terms(hpo, c("HP:0001873"), 
#'	list(term_sets=list("HP:0001873", "HP:0001902")), n=2)
#' @export
#' @import ontologyIndex
#' @import functional
#' @importFrom stats setNames
n_most_frequent_terms <- function(ontology, term_sets, n, terms=unique(unlist(term_sets))) {
	if (length(terms) <= n) {
		terms
	} else {
		counts <- Compose(unlist, table)(
			lapply(
				term_sets,
				get_ancestors,
				ontology=ontology
			)
		)

		inc.terms <- setNames(ifelse(terms %in% names(counts), counts[terms], 0), terms)

		cut.off <- sort(inc.terms, decreasing=TRUE)[n+1]
		names(inc.terms[inc.terms > cut.off])
	}
}

#' Remove terms with less than certain number of occurrences
#'
#' @template ontology
#' @template term_sets
#' @param n Integer
#' @template terms
#' @return Character vector
#' @seealso \code{\link{n_most_frequent_terms}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' remove_terms_with_less_than_n_occurrences(hpo, 
#'	term_sets=list("HP:0001873", "HP:0001902"), n=2)
#' @export
#' @import ontologyIndex
remove_terms_with_less_than_n_occurrences <- function(ontology, term_sets, n, terms=unique(unlist(term_sets))) intersect(
	terms,
	Compose(unlist, table, function(x) which(x >= n), names)(
		lapply(term_sets, get_ancestors, ontology=ontology)
	)
)

#' Remove terms which just link two other terms 
#'
#' @template ontology
#' @template terms
#' @param hard Logical value indicating whether to remove alternative direct paths to leaf nodes.
#' @return Character vector.
#' @seealso \code{\link{remove_uninformative_terms}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' remove_links(hpo, c("HP:0001873"))
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
remove_links <- function(ontology, terms, hard=FALSE) {
	repeat {
		adj <- get_pseudo_adjacency_matrix(ontology, terms)
		one.in.one.out <- apply(adj, 1, sum) <= 1 & apply(adj, 2, sum) == 1
		leaf <- apply(adj, 2, sum) == 0
		leaf.children.have.multiple.parents <- if (hard) FALSE else sapply(setNames(nm=terms), function(t) max(0, apply(adj[leaf & adj[,t],,drop=FALSE], 1, sum))) > 1
		to_remove <- one.in.one.out & !leaf.children.have.multiple.parents 
		if (sum(to_remove) > 0) {
			terms <- terms[!to_remove]
		} else {
			return(terms)
		}
	}
}

#' Get p-values for observing at least as many of each term as occur in \code{term_sets} given the population frequencies of the terms 
#'
#' @template ontology
#' @template term_sets
#' @param terms_freq Numeric vector of population frequencies of terms.
#' @return Numeric vector of log p-values named by correspnding term.
#' @seealso \code{\link{width_by_significance}}
#' @export
#' @import ontologyIndex
#' @importFrom stats pbinom
p_values_for_occurrence_of_term_in_group <- function(ontology, term_sets, terms_freq) pbinom(
	q=sapply(
		names(terms_freq),
		function(term) sum(
			sapply(
				lapply(
					term_sets, 
					function(x) get_ancestors(ontology, x)
				), 
				function(ancs) term %in% ancs
			)
		)
	)-1,
	size=length(term_sets),
	prob=terms_freq,
	lower.tail=FALSE,
	log.p=TRUE
)

#' Function to scale values between two given limits
#'
#' Could be useful to modify a vector of sizes to between, say 1 and 3, before passing to `onto_plot`.
#'
#' @param x Numeric vector
#' @param high Numeric value of largest size
#' @param low Numeric value of smallest size
#' @return Numeric vector
#' @examples
#' calibrate_sizes(c("HP:0000001"=10, "HP:0000006"=5), high=3, low=1)
#' @export
#' @import ontologyIndex
calibrate_sizes <- function(x, high, low) "+"(
	low,
	"*"(
		"/"(
			x-min(x),
			"if"(
				diff(range(x)) == 0,
				1,
				diff(range(x))
			)
		),
		high-low
	)
)

#' Remove terms not descending from phenotypic abnormality
#'
#' @template ontology
#' @template terms
#' @return Character vector.
#' @seealso \code{\link{remove_terms_with_less_than_n_occurrences}}, \code{\link{n_most_frequent_terms}}
#' @export
#' @import ontologyIndex
only_phenotype_abnormalities <- function(ontology, terms) {
	pa <- ontology$id[ontology$name == "Phenotypic abnormality"]
	Filter(
		x=terms,
		f=function(x) "%in%"(
			pa,
			ontology$ancestors[[x]]
		)
	)
}

#' Get ontology_plot object
#'
#' @template ontology
#' @template term_sets
#' @template frequencies
#' @template terms 
#' @param edge_attributes List of properties to set for arrows (note, these properties will be used for all arrow).
#' @param fillcolor Character vector of colours to fill nodes corresponding to \code{terms} with. Alternatively a function to set the colours of the nodes in the graph based on \code{term_sets}.
#' @param label Character vector of labels (or function to set them).
#' @param color Character vector of colours for borders of nodes representing \code{terms} (or function to set them).
#' @param width Numeric vector of widths for nodes (of function to set them).
#' @param fontsize Numeric vector of font sizes for the text to be placed in the nodes (or function to set them).
#' @param style Display style for nodes, defaults to \code{"filled"}.
#' @param fixedsize Character indicating whether nodes should be fixed size, "true", or adjusted to fit around the contained text, "false".
#' @param shape Character vector of shape names for nodes (or function to set them). Defaults to \code{"circle"}.
#' @param ... Other node attributes for dot format.
#' @return \code{ontology_plot} object.
#' @seealso \code{\link{write_dot}}
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' hpo_phenotypes <- c(
#' 	A=c("HP:0001382","HP:0004272","HP:0007917","HP:0004912","HP:0001596"),
#' 	B=c("HP:0001382","HP:0004272","HP:0002165","HP:0004800","HP:0004912"),
#' 	C=c("HP:0004800","HP:0001382","HP:0004912","HP:0007917","HP:0008743"),
#' 	D=c("HP:0001257","HP:0001382","HP:0007917","HP:0012623","HP:0002165"),
#' 	E=c("HP:0007917","HP:0004800","HP:0004272","HP:0001596","HP:0002165")
#' )
#' 
#' onto_plot(
#' 	ontology=hpo,
#' 	term_sets=hpo_phenotypes
#' )
#' @import Rgraphviz
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
onto_plot <- function(
	ontology,
	term_sets=NULL,
	frequencies=NULL,
	terms=remove_uninformative_terms(ontology, lapply(term_sets, get_ancestors, ontology=ontology)),
	edge_attributes=list(color="#000000", lty="solid"),
	fillcolor="powderblue",
	label=simple_labels,
	color="transparent",
	width=0.75,
	fontsize=15,
	style="filled",
	fixedsize="true",
	shape="circle",
	...
) {
	adj.mat <- t(get_pseudo_adjacency_matrix(
		ontology,
		terms
	))

	attrs <- lapply(
		c(list(style=style, fixedsize=fixedsize, fontsize=fontsize, shape=shape, width=width, fillcolor=fillcolor, label=label, color=color), list(...)),
		function(attribute) {
			switch(class(attribute),
				"function"=do.call(attribute, c(list(terms), Filter(f=Negate(is.null), lapply(setNames(nm=c("ontology", "term_sets", "frequencies")), function(argument.name) if (argument.name %in% names(formals(attribute))) get(argument.name))))),
				local({ right.length <- if (length(attribute) == 1) rep(attribute, times=length(terms)) else attribute; if (is.null(names(right.length))) setNames(right.length, terms) else right.length })
			) 
		}
	)

	structure(list(adjacency_matrix=adj.mat, node_attributes=attrs, edge_attributes=edge_attributes), class="ontology_plot")
}

#' Print function for \code{ontology_plot} object
#'
#' @param x Object of class ontologicalPlot.
#' @param ... Other options passed to be passed to plot().
#' @return Nothing. Side-effect: plots graphs.
#' @method print ontology_plot
#' @export
#' @import ontologyIndex
print.ontology_plot <- function(x, ...) {
	plot(x, ...)
}

#' Plotting function for \code{ontology_plot} object
#'
#' @param x Object of class ontologicalPlot.
#' @param ... Other options passed to plot().
#' @return Nothing, side-effect: plots a graph.
#' @method plot ontology_plot
#' @export
#' @import ontologyIndex
#' @importFrom methods new slot<-
plot.ontology_plot <- function(x, ...) {
	x$node_attributes$label <- gsub(x=x$node_attributes$label, pattern="\n", replacement="\\\\\n")
	hpo.graph <- new(
		"graphAM", 
		adjMat=x[["adjacency_matrix"]], 
		edgemode="directed"
	)

	result <- agopen(graph=hpo.graph, nodeAttrs=x[["node_attributes"]], name="ontological_plot") 
	if (length(result@AgEdge) > 0)
		for (i in 1:length(result@AgEdge)) {
			for (aai in 1:length(x[["edge_attributes"]])) slot(result@AgEdge[[i]], names(x[["edge_attributes"]])[aai]) <- x[["edge_attributes"]][[aai]]
		}

	plot(result, ...)
}

#' Function to write dot file
#'
#' @param ontology_plot Object of class `ontology_plot` to export.
#' @param file Character value of target file path.
#' @return Nothing, side effect - writes to file.
#' @seealso \code{\link{onto_plot}}
#' @export
#' @import ontologyIndex
write_dot <- function(ontology_plot, file) {
	n.terms <- nrow(ontology_plot$adjacency_matrix)
	node_lines <- paste(1:n.terms, " [", sapply(1:n.terms, function(x) paste(names(ontology_plot$node_attributes),"=\"",sapply(ontology_plot$node_attributes, "[", x), "\"", sep="", collapse=", ")), "];", sep="")
	edge_atts <- paste("[", paste(collapse=", ", sep="", names(ontology_plot$edge_attributes), "=\"", ontology_plot$edge_attributes, "\""), "]", sep="")
	edge_lines <- apply(which(arr.ind=TRUE, ontology_plot$adjacency_matrix), 1, function(x) paste(x[1], " -> ", x[2], " ", edge_atts, ";", sep=""))
	writeLines(text=c("digraph {", node_lines, edge_lines, "}"), con=file)
}

#' Split up node labels across lines so they fit in nodes better
#'
#' @template ontology
#' @template terms
#' @param official_names Logical value indicating whether to use the exact names from the ontology. Otherwise, shortened, capitalised names are used.
#' @return Character vector.
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_node_friendly_long_names(hpo, c("HP:0001873", "HP:0011877"))
#' @export
#' @import ontologyIndex
get_node_friendly_long_names <- function(ontology, terms, official_names=FALSE) {
	reorglabs <- if (official_names) ontology$name[terms] else sapply(
		gsub(
			"(Abnormality of (the )?)|(Abnormal)", 
			"", 
			ontology$name[terms]
		),
		simple_cap
	)

	reorglabs <- sapply(
		reorglabs, 
		function(x) {
			words <- strsplit(x, split=" |-")[[1]]
			if (length(words) == 1)
				return(words)
			
			lines <- list(words[1])
			for (word.no in 2:length(words))
				if (nchar(paste(c(words[word.no], lines[[length(lines)]]), collapse=" ")) > 17)
					lines <- c(lines, words[word.no])
				else
					lines[[length(lines)]] <- c(lines[[length(lines)]], words[word.no])

			desc.lines <- paste(
				lapply(lines, function(line) paste(line, collapse=" ")),

				collapse="\n"
			)

			desc.lines
		}
	)

	reorglabs
}
