## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(dev="svg", fig.width=7, fig.height=7)

## ------------------------------------------------------------------------
library(ontologyIndex)
library(ontologySimilarity)
library(SimReg)
data(hpo)
set.seed(1)
terms <- get_ancestors(hpo, c(hpo$id[match(c("Abnormality of thrombocytes","Hearing abnormality"), 
	hpo$name)], sample(hpo$id, size=50)))

## ------------------------------------------------------------------------
hearing_abnormality <- hpo$id[match("Hearing abnormality", hpo$name)]
genotypes <- c(rep(TRUE, 3), rep(FALSE, 97))

#give all subjects 5 random terms and add 'hearing abnormality' for those with y_i=TRUE
phenotypes <- lapply(genotypes, function(y_i) minimal_set(hpo, c(
if (y_i) hearing_abnormality else character(0), sample(terms, size=5))))

## ------------------------------------------------------------------------
samples <- sim_reg(ontology=hpo, x=phenotypes, y=genotypes)
print(summary(samples), ontology=hpo)

## ------------------------------------------------------------------------
lit_sims <- ifelse(grepl(x=hpo$name, ignore=TRUE, pattern="hearing"), 10, 1)
names(lit_sims) <- hpo$name

## ------------------------------------------------------------------------
thrombocytes <- hpo$id[match("Abnormality of thrombocytes", hpo$name)]
literature_phenotype <- c(hearing_abnormality, thrombocytes)
info <- get_term_info_content(hpo, phenotypes)

lit_sims_resnik <- apply(exp(get_term_set_to_term_sims(
	get_term_sim_mat(hpo, info, method="resnik"), 
	literature_phenotype)), 2, mean)

## ------------------------------------------------------------------------
with_prior_samples <- sim_reg(
	ontology=hpo,
	x=phenotypes,
	y=genotypes,
	lit_sims=lit_sims_resnik
)

print(summary(with_prior_samples), ontology=hpo)

## ----eval=FALSE----------------------------------------------------------
#  annotation_df <- read.table(header=FALSE, skip=1, sep="\t",
#  	file="ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt", stringsAsFactors=FALSE, quote="")
#  hpo_by_gene <- lapply(split(f=annotation_df[,2], x=annotation_df[,4]),
#  	function(trms) minimal_set(hpo, intersect(trms, hpo$id)))

