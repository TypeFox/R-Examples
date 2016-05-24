#' Function to set colours of nodes in plot to distinguish terms belonging to different term sets
#'
#' @template ontology
#' @template terms 
#' @template term_sets
#' @return Character vector of colours, named by term.
#' @seealso \code{\link{colour_by_frequency}}, \code{\link{colour_by_population_frequency}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
#' @importFrom grDevices rainbow
colour_by_term_set <- function(ontology, terms, term_sets) {
	ancestors.by.patient <- lapply(term_sets, get_ancestors, ontology=ontology)

	term.pat.mat <- data.frame(t(get_terms_by_set_matrix(ancestors.by.patient)))

	patient.combos <- unique(term.pat.mat)

	colours <- rainbow(nrow(patient.combos), alpha=0.5)
	
	setNames(
		colours[
			apply(term.pat.mat, 1, function(term.patients) which(apply(patient.combos, 1, function(combo) identical(term.patients, combo))))
		],
		terms
	)
}

#' Function to label nodes by \code{term_set}
#'
#' @template ontology
#' @template terms 
#' @template term_sets
#' @return Character vector of colours, named by term.
#' @seealso \code{\link{simple_labels}}, \code{\link{label_by_frequency}}, \code{\link{long_labels}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
label_by_term_set <- function(ontology, terms, term_sets) setNames(
	paste(
		get_node_friendly_long_names(ontology, terms),
		sapply(
			terms,
			function(term) paste(
				names(
					Filter(
						x=term_sets,
						f=function(patient.terms) term %in% get_ancestors(
							ontology,
							patient.terms
						)
					)
				),
				collapse=","
			)
		),
		sep="\n"
	),
	terms
)

#' Function to assign detailed node labels to terms
#'
#' Label includes term ID, term name, number of instances of term amongst \code{term_sets} and percentage frequency in population.
#'
#' @template ontology
#' @template terms
#' @template term_sets
#' @template frequencies
#' @return Character vector of labels, named by term.
#' @seealso \code{\link{simple_labels}}, \code{\link{label_by_frequency}}, \code{\link{label_by_term_set}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
long_labels <- function(ontology, terms, term_sets, frequencies) {

	node.friendly <- get_node_friendly_long_names(ontology, terms)
	freqs <- sapply(
		terms,
		function(term) sum(
			sapply(lapply(term_sets, function(x) get_ancestors(ontology, x)), function(ancs) term %in% ancs)
		)
	)

	result <- paste(
		terms,
		node.friendly,
		paste(
			round(100 * frequencies),	
			"% frequency", 
			sep=""
		),
		paste(
			freqs,
			" / ",
			length(term_sets),
			sep=""
		),
		sapply(
			terms,
			function(term) paste(
				names(
					Filter(
						x=term_sets,
						f=function(patient.terms) term %in% get_ancestors(
							ontology,
							patient.terms
						)
					)
				),
				collapse=","
			)
		),
		sep="\n"
	)
	
	setNames(
		result,
		terms
	)
}

#' Function to assign colours to terms based on frequency with which terms appear in \code{term_sets}
#'
#' @template ontology
#' @template terms 
#' @template term_sets
#' @param colour_func Function capable of returning a set of colours, given the number of colours it needs to return
#' @return Character vector of colours, named by term
#' @seealso \code{\link{colour_by_term_set}}, \code{\link{colour_by_population_frequency}}
#' @export
#' @import ontologyIndex
#' @importFrom grDevices colorRampPalette
colour_by_frequency <- function(
	ontology, 
	terms, 
	term_sets, 
	colour_func=colorRampPalette(c("Yellow", "Green", "#0099FF"))
) {
	ancestors.by.patient <- lapply(term_sets, function(x) get_ancestors(ontology, x))
	
	patients.with.term.count <- sapply(
		terms,
		function(term) sum(
			sapply(ancestors.by.patient, function(ancs) term %in% ancs)
		)
	)

	node.colours <- colour_func(1+diff(range(patients.with.term.count)))[
		patients.with.term.count-min(patients.with.term.count)+1
	]
	names(node.colours) <- terms
	
	node.colours <- node.colours[which(names(node.colours) %in% terms)]

	node.colours
}

#' Function to assign colours to terms based on population frequency of terms
#'
#' @template ontology
#' @template terms
#' @template frequencies
#' @param colour_palette Character vector of colours for the different information contents of the terms to be plotted, going from rare to common
#' @param max_colour_freq Numeric value in [0, 1] giving the maximum frequency (to which the dullest color will be assigned)
#' @param min_colour_freq Numeric value in [0, 1] giving the minimum frequency (to which the brightest color will be assigned)
#' @return Character vector of colours, named by term
#' @seealso \code{\link{colour_by_term_set}}, \code{\link{colour_by_frequency}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
colour_by_population_frequency <- function(
	ontology, 
	terms, 
	frequencies, 
	colour_palette=colorRampPalette(c("Yellow", "Green", "#0099FF"))(10),
	max_colour_freq=max(terms_freq),
	min_colour_freq=min(terms_freq)
) {
	terms_freq=frequencies[terms]
	freq.groups <- cut(terms_freq, seq(from=min_colour_freq, to=max_colour_freq, by=(max_colour_freq-min_colour_freq)/(length(colour_palette)-1)), include.lowest=TRUE)

	if (diff(range(terms_freq)) == 0)
		node.colors <- rep(colour_palette[1], length(terms_freq))
	else
		node.colors <- colour_palette[as.integer(freq.groups)]

	setNames(
		node.colors,
		terms	
	)
}

#' Function to get node sizes for terms based on statistical significance of seeing at least this number of each term in \code{term_sets}
#'
#' @template ontology
#' @template terms 
#' @template term_sets
#' @template frequencies
#' @return Character vector of sizes, named by term
#' @seealso \code{\link{width_by_frequency}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
width_by_significance <- function(ontology, terms, term_sets, frequencies) setNames(
	calibrate_sizes(
		-p_values_for_occurrence_of_term_in_group(
			ontology,
			term_sets,
			frequencies[terms]
		),
		3,
		1
	),
	terms
)

#' Function to get node sizes for terms based on frequency in \code{term_sets}
#'
#' @template ontology
#' @template terms
#' @template term_sets
#' @return Character vector of sizes, named by term
#' @seealso \code{\link{width_by_significance}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
width_by_frequency <- function(ontology, terms, term_sets) {
	group.freq <- sapply(
		terms,
		function(term) sum(
			sapply(
				lapply(
					term_sets, 
					function(x) get_ancestors(ontology, x)
				), 
				function(ancs) term %in% ancs
			)
		)
	)

	setNames(
		(
			function(x, high, low) "+"(
				low,
				"*"(
					"/"(
						x-min(x),
						diff(range(x))
					),
					high-low
				)
			)		
		)(group.freq, 3, 1),
		terms
	)
}
	
#' Function to get plot labels for terms based on frequency in \code{term_sets} 
#'
#' @template ontology
#' @template terms
#' @template term_sets
#' @return Character vector of labels, named by term.
#' @seealso \code{\link{simple_labels}}, \code{\link{long_labels}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
label_by_frequency <- function(ontology, terms, term_sets) {
	node.friendly <- get_node_friendly_long_names(ontology, terms)
	freqs <- sapply(
		terms,
		function(term) sum(
			sapply(lapply(term_sets, function(x) get_ancestors(ontology, x)), function(ancs) term %in% ancs)
		)
	)

	result <- paste(
		node.friendly,
		paste(
			freqs,
			" / ",
			length(term_sets),
			sep=""
		),
		sep="\n"
	)
	
	setNames(
		result,
		terms
	)
}
	
#' Get simplified labels for terms
#'
#' @template ontology
#' @template terms
#' @return Character vector of labels, named by term.
#' @seealso \code{\link{official_labels}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
simple_labels <- function(ontology, terms) setNames(
	get_node_friendly_long_names(ontology, terms),
	terms
)

#' Get official names for terms
#'
#' @template ontology
#' @template terms
#' @return Character vector of labels, named by term.
#' @seealso \code{\link{simple_labels}}
#' @export
#' @import ontologyIndex
#' @importFrom stats setNames
official_labels <- function(ontology, terms) setNames(
	get_node_friendly_long_names(ontology, terms, official_names=TRUE),
	terms
)

