## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(ontologyIndex))
suppressPackageStartupMessages(library(ontologySimilarity))
data(hpo)
set.seed(1)

## ------------------------------------------------------------------------
#random set of terms with ancestors
terms <- get_ancestors(hpo, sample(hpo$id, size=30))

#set information content of terms (as if each term occurs with frequency `1/n`)
information_content <- get_term_info_content(hpo, term_sets=as.list(terms))

#similarity of term pairs
tsm <- get_term_sim_mat(hpo, information_content)

## ------------------------------------------------------------------------
phenotypes <- replicate(simplify=FALSE, n=5, expr=minimal_set(hpo, sample(terms, size=8)))

## ------------------------------------------------------------------------
sim_mat <- get_sim_mat(tsm, phenotypes)
sim_mat

## ------------------------------------------------------------------------
get_sim(sim_mat, 1:3)

## ------------------------------------------------------------------------
get_sim_p(sim_mat, 1:3)

## ------------------------------------------------------------------------
get_sim_to_profile_p(tsm, phenotypes[[1]], phenotypes[[2]])

