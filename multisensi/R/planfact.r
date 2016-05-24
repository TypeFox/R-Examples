# Multisensi R package ; file planfact.r (last modified: 2013-07-03) 
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
planfact <- function(nb.niv,make.factor=TRUE)
#===========================================================================
{
    ## genere un plan factoriel complet dans l'ordre lexicographique
    ## ENTREES
    ## nb.niv: vecteur compose des nombres de niveaux de chaque facteur
    ## SORTIE
    ## un "dataframe" dont les colonnes sont des facteurs au sens Splus

    N <- prod(nb.niv)

    n <- length(nb.niv)

    plan <- matrix(NA,N,n)
    plan <- as.data.frame(plan)

    for(i in 1:(n-1))
    {
        plan[,i] <- rep( 1:nb.niv[i],  rep( prod(nb.niv[(i+1):n]), nb.niv[i] ))

        if(make.factor==TRUE){

            plan[,i] <- factor(plan[,i])

        }
    }

    plan[,n] <- rep( 1:nb.niv[n],  prod(nb.niv[1:(n-1)])   )

    if(make.factor==TRUE){plan[,n] <- factor(plan[,n])}

    if(make.factor==FALSE){plan <- plan-1}

    plan <- as.data.frame(plan)
    names(plan) <- paste("V",1:n,sep="")

    plan
}

