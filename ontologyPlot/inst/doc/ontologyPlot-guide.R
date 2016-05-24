## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(dev="svg")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(ontologyIndex))
suppressPackageStartupMessages(library(ontologyPlot))
data(hpo)

## ------------------------------------------------------------------------
onto_plot(hpo, terms=c("HP:0001873", "HP:0011877"))

## ------------------------------------------------------------------------
onto_plot(hpo, terms=get_ancestors(hpo, c("HP:0001873", "HP:0011877")))

## ------------------------------------------------------------------------
onto_plot(hpo, terms=remove_links(hpo, get_ancestors(hpo, c("HP:0001873", "HP:0011877"))))

## ------------------------------------------------------------------------
terms <- remove_links(hpo, get_ancestors(hpo, c("HP:0001873", "HP:0011877")))
onto_plot(hpo, terms=terms, label=terms, fillcolor=rainbow(length(terms)))

## ------------------------------------------------------------------------
onto_plot(hpo, terms=terms, label=official_labels)

## ------------------------------------------------------------------------
frequencies <- seq(from=0, to=1, by=1/length(terms))
names(frequencies) <- terms
onto_plot(hpo, terms=terms, frequencies=frequencies, 
	fillcolor=colour_by_population_frequency)

## ------------------------------------------------------------------------
hpo_phenotypes <- list(
	A=c("HP:0001382","HP:0004272","HP:0007917","HP:0004912","HP:0001596"),
	B=c("HP:0001382","HP:0004272","HP:0002165","HP:0004800","HP:0004912"),
	C=c("HP:0004800","HP:0001382","HP:0004912","HP:0007917","HP:0008743")
)

onto_plot(hpo, term_sets=hpo_phenotypes, label=label_by_term_set)

## ------------------------------------------------------------------------
onto_plot(
	hpo, 	
	frequencies=get_term_frequencies(hpo, hpo_phenotypes),
	term_sets=hpo_phenotypes,
	label=label_by_term_set,
	fillcolor=colour_by_frequency)

## ------------------------------------------------------------------------
onto_plot(
	hpo, 	
	frequencies=get_term_frequencies(hpo, hpo_phenotypes),
	term_sets=hpo_phenotypes,
	label=label_by_term_set,
	edge_attributes=list(color="red", lty="dashed"))

