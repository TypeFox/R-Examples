# Multisensi R package ; file simulmodel.r (last modified: 2016-01-07) 
# Copyright INRA 2011-2015 
# Authors: C. Bidot, M. Lamboni, H. Monod
# MaIAGE, INRA, Univ. Paris-Saclay, 78350 Jouy-en-Josas, France
#
# More about multisensi in http://cran.r-project.org/web/packages/multisensi/
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
#===========================================================================
simulmodel <- function(model,plan,nomFic=NULL,verbose=FALSE,...)
#===========================================================================
{
    ## simulations des sorties du modele selon le plan

    ##ENTREES
    ## nomFic est le nom du fichier ou se trouve le modele
    ##        (A entrer entre guillemets).
    ## model est la fonction renvoyant la sortie du
    ##         modele(Pas de guillemets).
    ## plan est le plan d'experience ou d'echantillonnage sous forme de data.frame


    ##SORTIES
    ##matrice de simulations des sortie sous forme array:  V
    ##Plan de simulations plan sous forme de data.frame : plan

    ## Chargement du fichier contenant le modele

    if(!is.null(nomFic)){ source(nomFic)}

    ## Renvoi de toutes les sorties du modele sous forme de liste.
    ## La taille de la liste est la meme
    ## que celle du plan factoriel, et chaque element de la liste
    ## est la sortie du modele (de la taille
    ## de la serie temporelle renvoyee par le modele).

    U <- NULL

    for (i in 1:nrow(plan)) {
        if( (verbose==TRUE) & ((i%%100)==0) ){ cat(i) }
        U <- c(U,list(model(plan[i,],...)))
    }
    if(verbose){cat("\n")}

    ## Transformation de U (format liste) en V sous forme de tableau pour l'ACP

    U2 <- unlist(U)
    V <- array(U2,dim=c(length(U[[1]]),nrow(plan)))
    V <- t(V)
    V <- as.data.frame(V)
    colnames(V)=names(U[[1]])
    return(V)
}

