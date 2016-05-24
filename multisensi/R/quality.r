# Multisensi R package ; file quality.r (last modified: 2015-09-14) 
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
quality <- function(echsimul,echsimul.app)
#===========================================================================
{
  ## calcule le biais, le coef. de determination, le lof

  ## INPUTS
  ## echsimul     : ARRAY des sorties ou matrice des sorties Y
  ## echsimul.app : sortie du modele approxime 

  ## OUTPUTS
  ## moy.biais    : biais integre (moyenne par colonne des biais)
  ## coef.det     : coefficient de determination R2
  ## residuals    : biais

  ## calcul de la moyenne de Y initiales
  moy.Y <- colMeans(echsimul)
  #colMeans et colSums + rapides que apply

  ## moyenne de Y sous forme matricielle
  moy.Ymat <- t(matrix(rep(moy.Y,nrow(echsimul)),ncol(echsimul),nrow(echsimul)))

  ## calcul du biais
  biais <- echsimul-echsimul.app
  moy.biais <- colMeans(biais)

  ##calcul du coefficient de determination
  ##  = 1 - SSresidus/SStotal
  coef.det <- 1 - colSums((biais)^2)/colSums((echsimul-moy.Ymat)^2)

  ## renvoi des sorties
  return(list(moy.biais=moy.biais,coef.det=coef.det,residuals=biais))
}

