
# correspond au fichier fonction_scores_covs_v2.R (fonction utilitaire a la fin de ce fichier) dans le dossier "programmes"

converti.terme = function(vec,n.loc)
{
il = as.numeric(vec)
ifelse (length(il)>1, n.loc + ifelse(il[1]>1,(il[1] - 1)*(n.loc - il[1]/2),0) + il[2] - il[1],il)     
}     
                 