## ----echo=FALSE----------------------------------------------------------
suppressPackageStartupMessages(library(gsEasy))
set.seed(1)

## ------------------------------------------------------------------------
#highly enriched... the set of ranks are relatively high out of 1000
gset(S=1:5 * 2, N=1000)

#random sets...
replicate(n=10, expr=gset(S=sample.int(n=1000, size=5), N=1000))

## ------------------------------------------------------------------------
gset(S=c("gene 1", "gene 5", "gene 40"), r=paste("gene", 1:100))

## ------------------------------------------------------------------------
gene_sets <- c(list(1:5), replicate(n=10, simplify=FALSE, expr=sample.int(n=1000, size=5)))
names(gene_sets) <- c("enriched set", paste("unenriched set", 1:10))
gene_sets
sapply(gene_sets, function(set) gset(S=set, N=1000))

## ------------------------------------------------------------------------
library(ontologyIndex)
data(hpo)
df <- data.frame(
	gene=c("gene 1", "gene 2"), 
	term=c("HP:0000598", "HP:0000118"), 
	name=hpo$name[c("HP:0000598", "HP:0000118")], 
	stringsAsFactors=FALSE,
	row.names=NULL)
df
get_ontological_gene_sets(hpo, gene=df$gene, term=df$term)

