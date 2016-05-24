# Multisensi R package ; file grpe.gsi.r (last modified: 2013-07-03) 
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
grpe.gsi <- function(GSI,fact.interet)
#===========================================================================
{
    ## GSI           : GSI object
    ## Fact.interet  ; a vecteur of the names of factors to be grouped.  Eg c("A","B")


    ##     fact.interet1 <- fact.interet
    ## Decomposition entre effets d'interet et autres
    facteurs <- rownames(GSI$Att)
    fact.interet <- apply(outer(facteurs,fact.interet,"=="),1,any)
    if(sum(fact.interet)==1){
        filtre1 <- GSI$Att[fact.interet,]
    } else{
        filtre1 <- apply( GSI$Att[fact.interet,],2,any)
    }
    if(sum(!fact.interet)==1){
        filtre2 <- GSI$Att[!fact.interet,]
    }  else{
        filtre2 <- apply( GSI$Att[!fact.interet,],2,any)
    }
    term1 <- filtre1 & (!filtre2)
    term2 <- (!filtre1) & filtre2
    indic.intrt <- rbind(GSI$Att[!fact.interet,],fact.interet=term1)
    ##     rownames(indic.intrt)[nrow( indic.intrt)] <- quote(fact.interet1)
    intrt.ss <- indic.intrt%*%as.matrix(GSI$SI)

    intrt.ss
}

