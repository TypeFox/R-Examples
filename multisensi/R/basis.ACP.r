# Multisensi R package ; file basis.ACP.r (last modified: 2016-02-03) 
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
basis.ACP <- function(simuls, basis.args=list())
#===========================================================================
{
  ## fait l'ACP, utilise prcomp()

  ## INPUTS
  ## simuls     : data.frame sur laquelle faire l'ACP (sortie des simulations du modele)
  ## basis.args : arguments specifiques a la methode de decomposition sous forme de list
  #               vide pour l'ACP

  ## OUTPUTS
  ## H          : coefficients/composantes principales de l'ACP, sortie x de la fonction prcomp
  ##              de taille nrow(simuls) x ncol(simuls)
  ## L          : matrice des vecteurs de la base ACP, sortie rotation de la fonction prcomp
  ##              de taille ncol(simuls) x ncol(simuls)
  ##  NB : simuls = H %*% t(L)
  ## call.info  : contient des indications sur la fonction utilisee, call.info$reduction="pca"


########## ATTENTION
# mettre en place des tests sur basis.args pour arguments optionnels de prcomp ?...
#
#ACP.prcomp <- do.call(prcomp,c(x=simuls,scale = normalized,basis.args))
#=> Erreur dans as.matrix(x) : 
#  l'argument "x" est manquant, avec aucune valeur par défaut
####################################################

  ACP.prcomp <- prcomp(simuls, center=FALSE, scale.=FALSE)
  # l'objet simuls doit etre centre-reduit si on veut faire ACP centree-reduite
  # pretraitement prevu dans fonction multivar

  return(list(H=ACP.prcomp$x,         # taille nrow(simuls) x ncol(simuls)
              L=ACP.prcomp$rotation,  # taille ncol(simuls) x ncol(simuls)
              call.info=list(reduction="pca")
              )
        )
}

