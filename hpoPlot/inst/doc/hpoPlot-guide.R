## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(hpoPlot))
data(hpo.terms)

## ------------------------------------------------------------------------
hpo.plot(hpo.terms, terms=c("HP:0001873", "HP:0011877"))

## ------------------------------------------------------------------------
hpo.plot(hpo.terms, terms=get.ancestors(hpo.terms, c("HP:0001873", "HP:0011877")))

## ------------------------------------------------------------------------
hpo.plot(hpo.terms, terms=remove.links(hpo.terms, get.ancestors(hpo.terms, c("HP:0001873", "HP:0011877"))))

## ------------------------------------------------------------------------
terms <- remove.links(hpo.terms, get.ancestors(hpo.terms, c("HP:0001873", "HP:0011877")))
hpo.plot(hpo.terms, terms=terms, colours=rainbow(length(terms)), labels=terms)

## ------------------------------------------------------------------------
hpo.plot(hpo.terms, terms=terms, labels=get.code.node.labels)

## ------------------------------------------------------------------------
plotting.context <- list(frequency=seq(from=0, to=1, by=1/length(terms)))
names(plotting.context$frequency) <- terms
hpo.plot(hpo.terms, terms=terms, plotting.context=plotting.context, 
	labels=get.code.node.labels, colours=get.pop.frequency.based.colours)

## ------------------------------------------------------------------------
phenotype.strings <- c(
	A="HP:0001382,HP:0004272,HP:0007917,HP:0004912,HP:0001596",
	B="HP:0001382,HP:0004272,HP:0002165,HP:0004800,HP:0004912",
	C="HP:0004800,HP:0001382,HP:0004912,HP:0007917,HP:0008743"
)

hpo.phenotypes <- term.set.list.from.character(phenotype.strings)

hpo.plot(hpo.terms, plotting.context=list(hpo.phenotypes=hpo.phenotypes), labels=get.case.based.labels)

## ------------------------------------------------------------------------
reduced.terms <- remove.uninformative.terms(hpo.terms, hpo.phenotypes)

hpo.plot(hpo.terms, plotting.context=list(hpo.phenotypes=hpo.phenotypes),
	terms=reduced.terms,
	labels=get.case.based.labels)

## ------------------------------------------------------------------------
plotting.context <- list(frequency=get.term.frequencies(hpo.terms, hpo.phenotypes),
	hpo.phenotypes=hpo.phenotypes)

hpo.plot(hpo.terms, plotting.context=list(hpo.phenotypes=hpo.phenotypes),
	terms=reduced.terms,
	labels=get.case.based.labels,
	colours=get.frequency.based.colours)

