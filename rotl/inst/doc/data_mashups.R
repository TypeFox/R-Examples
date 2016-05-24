## ---- data---------------------------------------------------------------
csv_path <- system.file("extdata", "protist_mutation_rates.csv", package = "rotl")
mu <- read.csv(csv_path, stringsAsFactors=FALSE)
mu

## ---- context------------------------------------------------------------
library(rotl)
tnrs_contexts()

## ---- match--------------------------------------------------------------
taxon_search <- tnrs_match_names(names=mu$species, context_name="All life")
knitr::kable(taxon_search)

## ---- munge--------------------------------------------------------------
mu$ott_name <- taxon_search$unique_name
mu$ott_id <- taxon_search$ott_id

## ---- properties---------------------------------------------------------
studies_properties()

## ----taxon_count---------------------------------------------------------
studies_find_trees(property="ot:ottId", value="180195")

## ---- all_taxa_count-----------------------------------------------------
hits <- lapply(mu$ott_id, studies_find_trees, property="ot:ottId", detailed = FALSE)
sapply(hits, function(x) sum(x[["n_matched_trees"]]))

## ----subtree,  fig.width=7, fig.height=4---------------------------------
tr <- tol_induced_subtree(ott_ids=mu$ott_id)
plot(tr)

## ---- match_names--------------------------------------------------------
mu$ott_name[1]
tr$tip.label[4]

## ---- sub----------------------------------------------------------------
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores=TRUE)
tr$tip.label %in% mu$ott_name

## ----phylobase-----------------------------------------------------------
library(phylobase)
mu_numeric <- mu[,c("mu", "pop.size", "genome.size")]
rownames(mu_numeric) <- mu$ott_name
tree_data <- phylo4d(tr, mu_numeric)

## ----  fig.width=7, fig.height=5-----------------------------------------
plot(tree_data)

## ------------------------------------------------------------------------
extra_data <- study_external_IDs("pg_1980")
extra_data

## ------------------------------------------------------------------------
library(rentrez)
seqs <- entrez_fetch(db="nucleotide", id=extra_data$nucleotide_ids[1:2], rettype="fasta")
cat(seqs)

## ------------------------------------------------------------------------
Tt_ids <- taxon_external_IDs(mu$ott_id[2])
Tt_ids

