#' Select most frequently annotated terms from a set of phenotypes
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @param n Integer 
#' @return Character vector of length at most \code{n}
#' @seealso \code{\link{remove.terms.with.less.than.n.occurrences}}, \code{\link{remove.uninformative.for.plot}}
#' @examples
#' data(hpo.terms)
#' n.most.frequent.terms(hpo.terms, c("HP:0001873"), 
#'	list(hpo.phenotypes=list("HP:0001873", "HP:0001902")), n=2)
#' @export
#' @import magrittr
#' @import functional
n.most.frequent.terms <- function(hpo.terms, terms, plotting.context, n) intersect(
	terms,
	Compose(	
		unlist,
		table, 
		function(x) sort(x, decreasing=TRUE), 
		names
	)(
		lapply(
			plotting.context$hpo.phenotypes,
			function(x) get.ancestors(hpo.terms, x)
		)
	)[1:min(n, length(terms))]
)

#' Remove terms with less than certain number of occurrences
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @param n Integer
#' @return Character vector
#' @seealso \code{\link{n.most.frequent.terms}}, \code{\link{remove.uninformative.for.plot}}
#' @examples
#' data(hpo.terms)
#' remove.terms.with.less.than.n.occurrences(hpo.terms, 
#'	c("HP:0001873"), list(hpo.phenotypes=list("HP:0001873", "HP:0001902")), 2)
#' @export
#' @import magrittr
remove.terms.with.less.than.n.occurrences <- function(hpo.terms, terms, plotting.context, n) intersect(
	terms,
	Compose(unlist, table, function(x) which(x >= n), names)(
		lapply(
			plotting.context$hpo.phenotypes,
			function(x) get.ancestors(hpo.terms, x)
		)
	)
)

#' Remove uninformative terms (fitting plotting filter format)
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector
#' @seealso \code{\link{remove.terms.with.less.than.n.occurrences}}, \code{\link{n.most.frequent.terms}}
#' @examples
#' data(hpo.terms)
#' remove.uninformative.for.plot(hpo.terms, 
#'	c("HP:0001873"), list(hpo.phenotypes=list("HP:0001873", "HP:0001902")))
#' @export
#' @import magrittr
remove.uninformative.for.plot <- function(hpo.terms, terms=get.ancestors(hpo.terms, unlist(plotting.context$hpo.phenotypes)), plotting.context) intersect(
	remove.uninformative.terms(
		hpo.terms,
		plotting.context$hpo.phenotypes
	),
	terms
)

#' Remove terms with exactly one parent and child from plot
#'
#' @template terms
#' @template hpo.terms
#' @template plotting.context
#' @return Character vector
#' @seealso \code{\link{remove.terms.with.less.than.n.occurrences}}, \code{\link{n.most.frequent.terms}}
#' @examples
#' data(hpo.terms)
#' remove.links(hpo.terms, c("HP:0001873"), list(hpo.phenotypes=list("HP:0001873", "HP:0001902")))
#' @export
#' @import magrittr
remove.links <- function(hpo.terms, terms, plotting.context=NULL) Filter(
	x=terms,
	f=function(term) !(
		length(intersect(terms, hpo.terms$children[[term]])) == 1
	)
)

#' Remove terms not descending from phenotypic abnormality
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector
#' @seealso \code{\link{remove.terms.with.less.than.n.occurrences}}, \code{\link{n.most.frequent.terms}}
#' @export
#' @import magrittr
remove.non.pa.terms <- function(hpo.terms, terms, plotting.context) Filter(
	x=terms,
	f=Curry(
		function(is.this, a.descendant.of.this) "%in%"(
			a.descendant.of.this,
			get.ancestors(
				hpo.terms,
				is.this	
			)
		),
		a.descendant.of.this=hpo.terms$id[hpo.terms$name == "Phenotypic abnormality"]
	)
)

#' Apply a list of term filters to a given plotting context
#'
#' @template hpo.terms
#' @template plotting.context
#' @param term.filters List of term filtering functions
#' @param starting.terms Character vector of HPO term codes to filter. Defaults to all terms in the `hpo.phenotypes' element of plotting.context, if it is present
#' @return Character vector of terms
#' @examples
#' data(hpo.terms)
#' apply.term.filters(hpo.terms=hpo.terms, plotting.context=list(
#'	hpo.phenotypes=list(Case1="HP:0001873")), term.filters=list(remove.links))
#' @export
#' @import magrittr
apply.term.filters <- function(hpo.terms, plotting.context, term.filters, starting.terms=NULL) Reduce(
	f=function(terms, filter.func) filter.func(terms=terms, hpo.terms=hpo.terms, plotting.context=plotting.context),
	init="if"(
		is.null(starting.terms),
		Compose(unlist, unique)(
			lapply(
				plotting.context$hpo.phenotypes,
				function(x) get.ancestors(hpo.terms, x)
			)
		),
		starting.terms
	),
	x=term.filters
)

#' Get HPO graph object
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @param colours Function to set the colours of the HPO nodes in the graph based on the plotting context, or a character vector of colours
#' @param labels Function to set the labels of the HPO nodes in the graph based on the plotting context, or a character vector of node labels
#' @param borders Function to set the borders of the HPO nodes in the graph based on the plotting context, or a character vector of border colours
#' @param sizes Function to set the sizes of the HPO nodes in the graph based on the plotting context, or a numeric vector of node sizes
#' @param font.sizes Function to set the font sizes of the text to be placed in the HPO nodes in the graph based on the plotting context, or an integer vector of font sizes
#' @param shapes Function to set the shapes of the HPO nodes in the graph based on the plotting context, or a character vector of shape names (defaults to 'circle')
#' @param nodeAttrs Pass nodeAttrs directly to rgraphviz plotting function
#' @param arrowAttrs List of properties to set for arrows (note, these properties will be used for all arrow)
#' @return graphAM S4 object 
#' @seealso \code{\link{hpo.plot}}
#' @examples
#' data(hpo.terms)
#' phenotype.strings <- c(
#' 	A="HP:0001382,HP:0004272,HP:0007917,HP:0004912,HP:0001596",
#' 	B="HP:0001382,HP:0004272,HP:0002165,HP:0004800,HP:0004912",
#' 	C="HP:0004800,HP:0001382,HP:0004912,HP:0007917,HP:0008743",
#' 	D="HP:0001257,HP:0001382,HP:0007917,HP:0012623,HP:0002165",
#' 	E="HP:0007917,HP:0004800,HP:0004272,HP:0001596,HP:0002165"
#' )
#' 
#' hpo.phenotypes <- term.set.list.from.character(phenotype.strings)
#' 
#' get.hpo.graph(
#' 	hpo.terms=hpo.terms,
#' 	plotting.context=list(hpo.phenotypes=hpo.phenotypes)
#' )
#' @import Rgraphviz
#' @export
#' @import magrittr
get.hpo.graph <- function(
	hpo.terms,
	terms=apply.term.filters(hpo.terms=hpo.terms, plotting.context=plotting.context, term.filters=list()),
	plotting.context=NULL,
	colours="white",
	labels=get.simple.node.labels,
	borders="#FFFFFF00",
	sizes=0.75,
	font.sizes=rep(30, length(terms)),
	shapes=rep("circle", length(terms)),
	nodeAttrs=NULL,
	arrowAttrs=list(color="#000000")
) {
	adj.mat <- get.term.pseudo.adjacency.matrix(
		hpo.terms,
		terms
	)

	hpo.graph <- new(
		"graphAM", 
		adjMat=adj.mat %>% t, 
		edgemode="directed"
	)

	if (is.null(nodeAttrs)) {
		attrs <- list(
			fontsize=if (class(font.sizes) == "function") font.sizes(hpo.terms, terms, plotting.context) else font.sizes %>% (function(x) if (length(x) == 1) rep(x, times=length(terms)) else x) %>% (function(x) if (is.null(names(x))) { names(x) <- terms; x } else x),
			shape=if (class(shapes) == "function") shapes(hpo.terms, terms, plotting.context) else shapes %>% (function(x) if (length(x) == 1) rep(x, times=length(terms)) else x) %>% (function(x) if (is.null(names(x))) { names(x) <- terms; x } else x),
			width=if (class(sizes) == "function") sizes(hpo.terms, terms, plotting.context) else sizes %>% (function(x) if (length(x) == 1) rep(x, times=length(terms)) else x) %>% (function(x) if (is.null(names(x))) { names(x) <- terms; x } else x),
			color=if (class(borders) == "function") borders(hpo.terms, terms, plotting.context) else borders %>% (function(x) if (length(x) == 1) rep(x, times=length(terms)) else x) %>% (function(x) if (is.null(names(x))) { names(x) <- terms; x } else x),
			fillcolor=if (class(colours) == "function") colours(hpo.terms, terms, plotting.context) else colours %>% (function(x) if (length(x) == 1) rep(x, times=length(terms)) else x) %>% (function(x) if (is.null(names(x))) { names(x) <- terms; x } else x),
			label=if (class(labels) == "function") labels(hpo.terms, terms, plotting.context) else labels %>% (function(x) if (length(x) == 1) rep(x, times=length(terms)) else x) %>% (function(x) if (is.null(names(x))) { names(x) <- terms; x } else x) 
		)
	} else {
		attrs <- nodeAttrs
	}
	result <- agopen(graph=hpo.graph, nodeAttrs=attrs, name="HPO.term.plot") 
	if (length(result@AgEdge) > 0)
		for (i in 1:length(result@AgEdge)) {
			result@AgEdge[[i]]@lty <- ifelse(
				result@AgEdge[[i]]@head %in% hpo.terms$parents[[result@AgEdge[[i]]@tail]],
				"solid",
				"dashed"
			)
			for (aai in 1:length(arrowAttrs)) slot(result@AgEdge[[i]], names(arrowAttrs)[aai]) <- arrowAttrs[[aai]]
			
		}

	result
}

#' Plot HPO graph object
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @template hpo.phenotypes 
#' @param term.frequencies Numeric vector of population frequencies of terms (named by term codes)
#' @param colours Function to set the colours of the HPO nodes in the graph based on the plotting context, or a character vector of colours
#' @param labels Function to set the labels of the HPO nodes in the graph based on the plotting context, or a character vector of node labels
#' @param borders Function to set the borders of the HPO nodes in the graph based on the plotting context, or a character vector of border colours
#' @param sizes Function to set the sizes of the HPO nodes in the graph based on the plotting context, or a numeric vector of node sizes
#' @param font.sizes Function to set the font sizes of the text to be placed in the HPO nodes in the graph based on the plotting context, or an integer vector of font sizes
#' @param shapes Function to set the shapes of the HPO nodes in the graph based on the plotting context, or a character vector of shape names (defaults to 'circle')
#' @param nodeAttrs Pass nodeAttrs directly to rgraphviz plotting function
#' @param arrowAttrs List of properties to set for arrows (note, these properties will be used for all arrow)
#' @param ... Extra arguments to pass to plot
#' @return Plots graph
#' @seealso \code{\link{get.hpo.graph}}
#' @examples
#' data(hpo.terms)
#' hpo.plot(
#' 	hpo.terms=hpo.terms,
#' 	terms=get.ancestors(hpo.terms, 
#'		c("HP:0001382","HP:0004272","HP:0007917","HP:0004912","HP:0001596"))
#' )
#' @export
#' @import magrittr
hpo.plot <- function(
	hpo.terms,
	terms=apply.term.filters(hpo.terms=hpo.terms, plotting.context=plotting.context, term.filters=list()),
	plotting.context=NULL,
	hpo.phenotypes=NULL,
	term.frequencies=NULL,
	colours="cyan",
	labels=get.simple.node.labels,
	borders="#FFFFFF00",
	sizes=0.75,
	font.sizes=rep(30, length(terms)),
	shapes=rep("circle", length(terms)),
	nodeAttrs=NULL,
	arrowAttrs=list(color="#000000"),
	...
) {
	original.mar <- par("mar")

	if (!is.null(hpo.phenotypes))
		plotting.context$hpo.phenotypes <- hpo.phenotypes
	if (!is.null(term.frequencies))
		plotting.context$information <- -log(term.frequencies)

	graph <- get.hpo.graph(
		hpo.terms=hpo.terms,
		terms=terms,
		plotting.context=plotting.context,
		colours=colours,
		font.sizes=font.sizes,
		shapes=shapes,
		labels=labels,
		borders=borders,
		sizes=sizes,
		nodeAttrs=nodeAttrs,
		arrowAttrs=arrowAttrs
	)
	
	plot(
		graph,
		...
	)

	par(mar=original.mar)
}
	
#' Function to set colours of HPO nodes in plot to distinguish terms belonging to different sets of phenotypes
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @return Character vector of colours, named by term
#' @export
#' @import magrittr
get.case.based.colours <- function(hpo.terms, terms, plotting.context) {
	hpo.phenotypes <- plotting.context$hpo.phenotypes
	ancestors.by.patient <- lapply(hpo.phenotypes, function(x) get.ancestors(hpo.terms, x))

	term.pat.mat <- data.frame(t(get.case.term.matrix(ancestors.by.patient)))

	patient.combos <- unique(term.pat.mat)

	colours <- rainbow(nrow(patient.combos), alpha=0.5)
	
	setNames(
		colours[
			apply(term.pat.mat, 1, function(term.patients) which(apply(patient.combos, 1, function(combo) identical(term.patients, combo))))
		],
		terms
	)
}

#' Function to label HPO nodes in plot to indicate to which phenotypes each of the terms belong
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @return Character vector of colours, named by term
#' @export
#' @import magrittr
get.case.based.labels <- function(hpo.terms, terms, plotting.context) setNames(
	paste(
		terms,
		get.node.friendly.long.names(hpo.terms, terms),
		sapply(
			terms,
			function(term) paste(
				names(
					Filter(
						x=plotting.context$hpo.phenotypes,
						f=function(patient.terms) term %in% get.ancestors(
							hpo.terms,
							patient.terms
						)
					)
				),
				collapse=","
			)
		),
		sep="\\\n"
	),
	terms
)

#' Function to label HPO nodes in plot with full labels
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector of labels, named by term
#' @export
#' @import magrittr
get.full.labels <- function(hpo.terms, terms, plotting.context) {

	node.friendly <- get.node.friendly.long.names(hpo.terms, terms)
	freqs <- sapply(
		terms,
		function(term) sum(
			sapply(lapply(plotting.context$hpo.phenotypes, function(x) get.ancestors(hpo.terms, x)), function(ancs) term %in% ancs)
		)
	)

	result <- paste(
		terms,
		node.friendly,
		paste(
			round(100 * exp(-plotting.context$information[terms])),	
			"% Of Collection", 
			sep=""
		),
		paste(
			freqs,
			" / ",
			length(plotting.context$hpo.phenotypes),
			sep=""
		),
		sapply(
			terms,
			function(term) paste(
				names(
					Filter(
						x=plotting.context$hpo.phenotypes,
						f=function(patient.terms) term %in% get.ancestors(
							hpo.terms,
							patient.terms
						)
					)
				),
				collapse=","
			)
		),
		sep="\\\n"
	)
	
	setNames(
		result,
		terms
	)
}

#' Function to colour HPO nodes in plot with colours based on frequency with which terms appear in phenotypes
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @param colour.func Function capable of returning a set of colours, given the number of colours it needs to return
#' @return Character vector of colours, named by term
#' @export
#' @import magrittr
get.frequency.based.colours <- function(
	hpo.terms, 
	terms, 
	plotting.context, 
	colour.func=NULL
) {
	if (is.null(colour.func))
		colour.func <- colorRampPalette(c("Yellow", "Green", "#0099FF"))

	ancestors.by.patient <- lapply(plotting.context$hpo.phenotypes, function(x) get.ancestors(hpo.terms, x))
	
	patients.with.term.count <- sapply(
		terms,
		function(term) sum(
			sapply(ancestors.by.patient, function(ancs) term %in% ancs)
		)
	)

	node.colours <- colour.func(1+diff(range(patients.with.term.count)))[
		patients.with.term.count-min(patients.with.term.count)+1
	]
	names(node.colours) <- terms
	
	#strip out the terms not relevant to this plot...
	node.colours <- node.colours[which(names(node.colours) %in% terms)]

	node.colours
}

#' Function to colour HPO nodes in plot with colours based on information content/frequency of terms with respect to population
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @param colourPalette Character vector of colours for the different information contents of the terms to be plotted, going from rare to common
#' @param terms.freq Numeric vector of frequencies of terms in plot, named by term
#' @param max.colour.freq Numeric value in [0, 1] giving the maximum frequency (to which the dullest color will be assigned)
#' @param min.colour.freq Numeric value in [0, 1] giving the minimum frequency (to which the brightest color will be assigned)
#' @return Character vector of colours, named by term
#' @export
#' @import magrittr
get.pop.frequency.based.colours <- function(
	hpo.terms, 
	terms, 
	plotting.context, 
	colourPalette=colorRampPalette(c("Yellow", "Green", "#0099FF"))(10),
	terms.freq=if (is.null(plotting.context$frequency)) exp(-plotting.context$information[terms]) else plotting.context$frequency,
	max.colour.freq=max(terms.freq),
	min.colour.freq=min(terms.freq)
) {
	freq.groups <- cut(terms.freq, seq(from=min.colour.freq, to=max.colour.freq, by=(max.colour.freq-min.colour.freq)/(length(colourPalette)-1)), include.lowest=TRUE)

	if (diff(range(terms.freq)) == 0)
		node.colors <- rep(colourPalette[1], length(terms.freq))
	else
		node.colors <- colourPalette[as.integer(freq.groups)]

	setNames(
		node.colors,
		terms	
	)
}

#under the null that the distribution of the (binomial) occurrence (i.e. the 'success rate' in gcse speak) of each term in the group is the same as that of the population, calculates a p.value for each of the terms for which information content is supplied - i.e. a vector of information contents named by term - the paramter 'terms'...
#' Get p-values for observing at least as many of each term as have been in phenotypes given information content
#'
#' @template hpo.terms
#' @template hpo.phenotypes
#' @param terms.freq Numeric vector of population frequencies of terms
#' @return Numeric vector of log p-values named by correspnding term
#' @export
#' @import magrittr
p.values.for.occurrence.of.term.in.group <- function(hpo.terms, hpo.phenotypes, terms.freq) pbinom(
	q=sapply(
		names(terms.freq),
		function(term) sum(
			sapply(
				lapply(
					hpo.phenotypes, 
					function(x) get.ancestors(hpo.terms, x)
				), 
				function(ancs) term %in% ancs
			)
		)
	)-1,
	size=length(hpo.phenotypes),
	prob=terms.freq,
	lower.tail=FALSE,
	log.p=TRUE
)

#' Function to scale sizes of terms between two given limits
#'
#' @param x Numeric vector of term relative sizes named by term
#' @param high Numeric vector of largest size
#' @param low Numeric vector of smallest size
#' @return Numeric vector
#' @examples
#' calibrate.sizes(c("HP:0000001"=10, "HP:0000006"=5), high=3, low=1)
#' @export
#' @import magrittr
calibrate.sizes <- function(x, high, low) "+"(
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

#' Function to size HPO nodes in plot with colours based on significance of seeing this many of each term in phenotypes 
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @return Character vector of sizes, named by term
#' @export
#' @import magrittr
get.significance.based.sizes <- function(hpo.terms, terms, plotting.context) setNames(
	calibrate.sizes(
		-p.values.for.occurrence.of.term.in.group(
			hpo.terms,
			plotting.context$hpo.phenotypes,
			exp(-plotting.context$information[terms])
		),
		3,
		1
	),
	terms
)

#' Function to size HPO nodes in plot based on frequency of occurrence in phenotypes
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector of sizes, named by term
#' @export
#' @import magrittr
get.frequency.based.sizes <- function(hpo.terms, terms, plotting.context) {
	group.freq <- sapply(
		terms,
		function(term) sum(
			sapply(
				lapply(
					plotting.context$hpo.phenotypes, 
					function(x) get.ancestors(hpo.terms, x)
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
	
#' Function to label HPO nodes in plot based on frequency of occurrence in phenotypes
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector of labels, named by term
#' @export
#' @import magrittr
get.frequency.based.labels <- function(hpo.terms, terms, plotting.context) {
	node.friendly <- get.node.friendly.long.names(hpo.terms, terms)
	freqs <- sapply(
		terms,
		function(term) sum(
			sapply(lapply(plotting.context$hpo.phenotypes, function(x) get.ancestors(hpo.terms, x)), function(ancs) term %in% ancs)
		)
	)

	result <- paste(
		terms,
		node.friendly,
		paste(
			round(100 * exp(-plotting.context$information[terms])),	
			"% Of Cases", 
			sep=""
		),
		paste(
			freqs,
			" / ",
			length(plotting.context$hpo.phenotypes),
			sep=""
		),
		sep="\\\n"
	)
	
	setNames(
		result,
		terms
	)
}
	
#' Split up the HPO term descriptions so they fit in nodes for plot
#'
#' @template hpo.terms
#' @template terms
#' @return Character vector
#' @examples
#' data(hpo.terms)
#' get.node.friendly.long.names(hpo.terms, c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.node.friendly.long.names <- function(hpo.terms, terms) {
	#get rid of the "Abnormality of"s because they take up too much room...
	simpleCap <- function(x) {
		s <- strsplit(x, " ")[[1]]
		paste(
			toupper(substring(s, 1,1)), 
			substring(s, 2),
			sep="", collapse=" "
		)
	}

	reorglabs <- sapply(
		gsub(
			"(Abnormality of (the )?)|(Abnormal)", 
			"", 
			hpo.terms$name[terms]
		),
		simpleCap
	)

	#hack... if there are more than 4 words, then put the first 4 words on diff line...
	reorglabs <- sapply(
		reorglabs, 
		function(x) {
			words <- strsplit(x, split=" |-")[[1]]
			if (length(words) == 1)
				return(words)
			
			lines <- list(words[1])
			for (word.no in 2:length(words))
				#if adding the next word to the current line makes it too long, go to new line
				if (nchar(paste(c(words[word.no], lines[[length(lines)]]), collapse=" ")) > 17)
					lines <- c(lines, words[word.no])
				else
					lines[[length(lines)]] <- c(lines[[length(lines)]], words[word.no])

			desc.lines <- paste(
				lapply(lines, function(line) paste(line, collapse=" ")),
				collapse="\\\n"
			)

			desc.lines
		}
	)

	reorglabs
}

#' Function to label HPO nodes in plot with just node description
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector of labels, named by term
#' @export
#' @import magrittr
get.simple.node.labels <- function(hpo.terms, terms, plotting.context) setNames(
	get.node.friendly.long.names(hpo.terms, terms),
	terms
)

#' Function to label HPO nodes in plot with node description and information content
#'
#' @template hpo.terms
#' @template terms 
#' @template plotting.context
#' @return Character vector of labels, named by term
#' @export
#' @import magrittr
get.informative.node.labels <- function(hpo.terms, terms, plotting.context) setNames(
	paste(
		terms,
		get.node.friendly.long.names(hpo.terms, terms),
		paste(
			round(100 * exp(-plotting.context$information[terms])),	
			"%", 
			sep=""
		),
		sep="\\\n"
	), 
	terms
)

#' Function to label HPO nodes in plot with just HPO code
#'
#' @template hpo.terms
#' @template terms
#' @template plotting.context
#' @return Character vector of labels, named by term
#' @export
#' @import magrittr
get.code.node.labels <- function(hpo.terms, terms, plotting.context) setNames(terms, terms)

