# Multisensi R package ; file planfact.as.r (last modified: 2013-07-03) 
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
planfact.as <- function(input)
#===========================================================================
{
    ##genere le plan factoriel complet en utilisant les vraies
    ##       modalites des parametres

    ##ENTREES
    ## input      est une liste d'elements contenant les modalites
    ##            des differents parametres


    ##SORTIES
    ## plan factoriel complet sous forme de data frame contenant
    ## les vraies modalites des inputs

    nb.niv <-1:length(input)
    for (k in 1:length(nb.niv)) {
        nb.niv[k] <- length(input[[k]])
    }
    ## Creation du plan factoriel
    plan <- planfact(nb.niv,make.factor=FALSE)
    ## Recodage des facteurs avec les valeurs transmises dans input

    echSimul <- matrix(NA,nrow(plan),ncol(plan))

    for (j in 1:ncol(plan)) {
        tmp <- input[[j]] ## modalites du facteur j
        plj <- plan[,j]+1   ## niveaux du facteur j dans le plan
        tmp <- tmp[plj]
        echSimul[,j] <- tmp

    }


    ## Recuperation des noms des variables

    if( is.null(names(input))) {  colnames(echSimul) <- colnames(plan)}
    else {                     colnames(echSimul) <- names(input[1:ncol(plan)])}

    ## Transformation du plan en dataframe
    echSimul <- as.data.frame(echSimul)

    ## Renvoi de l'echantillon
    return(echSimul)
}

