# Multisensi R package ; file basis.bsplines.r (last modified: 2016-02-03) 
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
basis.bsplines <- function(simuls, basis.args=list(knots=5, mdegree=3))
#===========================================================================
{
  ## fait la decomposition sur une base B-splines definie par les parametres knots et mdegree
  ## utilise sesBsplinesNORM()

  ## INPUTS
  ## simuls     : data.frame sur laquelle faire la decomposition (sortie des simulations du modele)
  ## basis.args : arguments specifiques a la methode de decomposition sous forme de list
  ## - knots    : nombre de noeuds ou vecteur avec leurs positions (5 par defaut)
  ## - mdegree  : degre des Bsplines (3 par defaut)
  ## - x.coord  : coordonnees ou calculer les valeurs des Bsplines (optionnel, 1:ncol(simuls) par defaut)
  ##              vecteur de longueur ncol(simuls)

  ## OUTPUTS
  ## H          : coefficients associes a la base de decomposition
  ##              de taille nrow(simuls) x nbsplines (nbsplines = knots+mdegree)
  ## L          : matrice des vecteurs de la base B-splines
  ##              de taille ncol(simuls) x nbsplines
  ##  NB : simuls ~ H %*% t(L)
  ## call.info  : contient des indications sur la fonction utilisee, call.info$reduction="b-splines"



  b.args=basis.args
  # definition des valeurs par defaut des arguments
  if(is.null(basis.args$knots)){
    b.args$knots=5;
#    cat("Warning : basis.bsplines argument 'knots' not given,\n          default value 'knots=5' used.\n")
  }
  if(is.null(basis.args$mdegree)){
    b.args$mdegree=3;
#    cat("Warning : basis.bsplines argument 'mdegree' not given,\n          default value 'mdegree=3' used.\n")
  }
  if(is.null(basis.args$x.coord)){
    b.args$x.coord=seq(1,ncol(simuls));
#    cat("Warning : basis.bsplines argument 'x.coord' not given,\n          default value 'x.coord=1:",ncol(simuls),"' used.\n")
  }else if(length(b.args$x.coord)!=ncol(simuls)){
    stop("basis.bsplines argument 'x.coord' length not adequate,\n         length(x.coord) must be ",ncol(simuls)," (",length(b.args$x.coord)," used).\n")
  }

  ##calcul des splines normees
  spl=sesBsplinesNORM(x=b.args$x.coord,knots=b.args$knots,m=b.args$mdegree)

  # nb de splines = knots+mdegree
  # spl$projecteur de taille nbsplines x ncol(simuls) equivaut bsplines^(-1)
  # spl$bsplines de taille ncol(simuls) x nbsplines

  ## calcul de H (H=Y.L^-1):
  # t(L) %*% t(spl$projecteur) = I et spl$projecteur %*% L = I
  # H*t(L) = Y => H*t(L) * t(spl$projecteur) = Y * t(spl$projecteur) => H =  Y * t(spl$projecteur)
  coefspl2 = as.matrix(simuls)%*%t(spl$projecteur);
  # coefspl2 de taille nrow(simuls) x nbsplines

  return(list(H=coefspl2,     # taille nrow(simuls) x nbsplines
              L=spl$bsplines, # taille ncol(simuls) x nbsplines
              call.info=list(reduction="b-splines")
              )
        )
}

