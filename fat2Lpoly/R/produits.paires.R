
# correspond au fichier fonction_scores_covs_v2.R (fonction utilitaire a la fin de ce fichier) dans le dossier "programmes"

produits.paires <- function(vec) 
{
mat <- (vec%o%vec)
mat[upper.tri(mat)]
}
