
# Jordie Croteau
# 17 aout 2012
# fonction a laquelle on fournit les noms des fichiers de format merlin et les autres parametres pour executer l'analyse "extension 2-locus de GDT"
# n.levels peut avoir la valeur 2 ou 4 (si n.levels=2, il faut mettre design.constraint=design.dichotomique)

# correspond au fichier fonction_globale_donnees_SPAP_v4.R dans le dossier "programmes"


# modifie le 22 aout (par Jordie)

# modifie le 23 aout (par Jordie) pour mettre la partie de la boucle principale dans une fonction separee: fat2Lpoly.withinR.

# petites modifs le 27 aout.

# 4 avril 2013: valeur par defaut des arguments ibdfilenames et ibd.loci fixee a NULL.
#               Si un de ces 2 arguments prend la valeur NULL, alors la fonction fat2Lpoly.withinR fait le calcul des coefficients de kinship (a priori) au lieu des IBD.

fat2Lpoly=function(pedfilenames,datfilenames,freq.data,ibdfilenames=NULL,snp.names.mat,ibd.loci=NULL,joint.tests=NULL,contingency.file=FALSE,design.constraint,par.constrained,constraints,pairweights=calcule.poids.alphafixe,lc=NULL,alpha=NULL)
{
###################### Definition des arguments #####################################################################################
# pedfilenames : vecteur des noms de fichiers ped (un fichier par locus).  Les sujets inclus peuvent etre un sous-ensemble de ceux 
#                inclus dans les fichiers d'IDB.
# datfilenames : vecteur des noms de fichiers dat (un fichier par locus). 
# freq.data: vecteur des noms de fichiers freq (un fichier par locus). 
# 
# Tous ces fichiers doivent etre en format Merlin (voir http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html pour la description detaillee de ce format)
#
# ibdfilenames : vecteur des noms de fichiers d'IBD (un fichier par locus). Ces fichiers doivent aussi avoir le format des resultats d'IBD de Merlin.
# n.levels : nombre de categories du statut (le statut est la combinaison des 6e colonne (endophenotype) et 7e colonne (phenotype) des fichiers ped)
# snp.names.mat : matrice de une ou deux colonnes donnant les noms des SNPs (si une colonne) ou des paires de SNPs (si deux colonnes) a analyser
#                 (pour le moment, on suppose qu'il y a maximum 2 locus)
# ibd.loci : matrice des memes dimensions que snp.names.mat, donnant les noms respectifs des microsats (ou les positions) les plus pres des SNP,
#                      pour extraction des resultats d'IBD.
# joint.tests : list of vectors of numbers between 1 and the total number of parameters in 'design.constraint'. Each vector gives parameter indices to test the corresponding parameters jointly.
# contingency.file: if 'TRUE' (default is 'FALSE'), then a file called descriptive_statistics.txt is created and contingency tables with the numbers of subjects per level are progressively added to this file.
# design.constraint : fonction construisant les matrices de design a l'INTeRIEUR chaque categorie pour des contraintes specifiques a chaque categorie.
#                     Les matrices de design comprenant seulement les effets principaux des locus qui sont utilise pour le calcul des covariances
#                     sont aussi produite par cette fonction
# par.constrained : vecteur des indices du parametre implique dans chaque contrainte ENTRE des categories(longueur nc). Utilise dans design.polytomous.
# constraints: matrice (K-1) x nc specifiant des contraintes entre les parametres  ENTRE les modeles logistiques
#              pour differentes categories de la variable reponse, une contrainte par colonne. Utilise dans design.polytomous.
# lc: numerical identifier of the SNP (locus) on which to condition when testing model terms. Defaults to NULL, or no conditioning.
#####################################################################################################################################

if(is.null(ibd.loci)|is.null(ibdfilenames)) cat("\n","Warning: Either the argument ibd.loci or ibdfilenames was not specified. The kinship coefficients (multiplied by 2) will be used in the computation of the score statistics (instead of the expectation of the IBD probabilities).","\n")

# lecture des donnees de format merlin
ped.x.all=read.merlin.files(pedfilenames,datfilenames,freq.data,ibdfilenames)

# execution des tests pour les SNPs ou paires de SNPs
tests.loop=fat2Lpoly.withinR(ped.x.all,snp.names.mat,ibd.loci,contingency.file,design.constraint,par.constrained,constraints,pairweights=pairweights,lc=lc,alpha=alpha)

# calcul des scores et valeur-p des differents tests pour tous les SNPs testes
p.values.scores=get.scores.pvalues(tests.loop,joint.tests)

return(list(scores.covs.all.SNPs=tests.loop$scores.covs.all.SNPs,p.values.scores=p.values.scores,MA.table=ped.x.all$MA.table,y1=ped.x.all$y1.name,y2=ped.x.all$y2.name))
}
