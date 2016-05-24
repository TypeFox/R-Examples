# Multisensi R package ; file yapprox.r (last modified: 2015-10-15) 
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
yapprox <- function(multivar.obj, nbcomp=2, aov.obj)
#===========================================================================
{
  ## calcule les sorties du modele approxime

  ## INPUTS
  ## multivar.obj : "objet" multivar, sortie de la fonction multivar (liste)
  ## nbcomp       : nombre de composantes a utiliser pour le metamodele
  ## aov.obj      : objet/liste issu de la sortie de la fonction analysis.anoasg

  ## OUTPUTS
  ## echsimul.app : sorties Y approximees sous forme de matrice

  ##recuperation des nbcomp premieres lignes de la matrice inversee des vecteurs propres.
  inv.comp <-  t(multivar.obj$L)[1:nbcomp,,drop=FALSE]
  # multivar.obj$L de taille ncol(Y) x nbcomp++

  ##calcul des sorties approximees des variables (PC.predict)
  echsimul.app <- aov.obj$Hpredict %*% inv.comp
  # aov.obj$Hpredict de taille nrow(Y)*nbcomp

  ## reconstitution des  valeurs de Yapp du fait que les Y utilises pour multivar.obj etaient (ou non) reduits
  echsimul.app <- echsimul.app %*% diag(multivar.obj$scaling,ncol(echsimul.app),ncol(echsimul.app));
  # si normalise, alors scaling=sdY, sinon, scaling=1

  ## decentre
  echsimul.app <- multivar.obj$centering + echsimul.app
  colnames(echsimul.app)=rownames(multivar.obj$L)


  return(echsimul.app)
}

