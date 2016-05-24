# Multisensi R package ; file basis.poly.r (last modified: 2016-02-03) 
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
basis.poly <- function(simuls, basis.args=list(degree=3))
#===========================================================================
{
  ## fait la decomposition sur une base de polynomes definie par le parametre degree
  ## utilise poly()

  ## INPUTS
  ## simuls     : data.frame sur laquelle faire la decomposition (sortie des simulations du modele)
  ## basis.args : arguments specifiques a la methode de decomposition sous forme de list
  ## - degree   : degre maximal des polynomes (3 par defaut)
  ## - x.coord  : coordonnees ou calculer les valeurs des polynomes (optionnel, 1:ncol(simuls) par defaut)

  ## OUTPUTS
  ## H          : coefficients associes a la base de decomposition
  ##              de taille nrow(simuls) x nbpoly (nbpoly = degree+1)
  ## L          : matrice des vecteurs de la base de polynomes
  ##              de taille ncol(simuls) x nbpoly
  ##  NB : simuls ~ H %*% t(L)
  ## call.info  : contient des indications sur la fonction utilisee, call.info$reduction="polynomial"

  b.args=basis.args
  # definition des valeurs par defaut des arguments
  if(is.null(basis.args$degree)){
    b.args$degree=3;
#    cat("Warning : basis.poly argument 'degree' not given,\n          default value 'degree=3' used.\n")
  }
  if(is.null(basis.args$x.coord)){
    b.args$x.coord=seq(1,ncol(simuls));
#    cat("Warning : basis.poly argument 'x.coord' not given,\n          default value 'x.coord=1:",ncol(simuls),"' used.\n")
  }else if(length(b.args$x.coord)!=ncol(simuls)){
    stop("basis.poly argument 'x.coord' length not adequate,\n         length(x.coord) must be ",ncol(simuls)," (",length(b.args$x.coord)," used).\n")
  }


  ##calcul des polynomes
  pbase=poly(b.args$x.coord,degree=b.args$degree)
  # ces polynomes sont orthogonaux au poly constant mais il faut l'ajouter 
  # => 1/sqrt(N) pour etre norme 
  # cad colSums(pbase^2) = 1 1 ... 1
  pbase=cbind(rep(1/sqrt(nrow(pbase)),nrow(pbase)),pbase)
  colnames(pbase)[1]="0"
  
  colnames(pbase) <- paste("Deg",colnames(pbase),sep="")

  # pbase de taille ncol(simuls) x (degree +1)

  #projecteur = solve( t(pbase) %*% pbase) %*% t(pbase)
  # "solve" solves the equation ‘a %*% x = b’ for ‘x’ 
  # ici a =t(pbase) %*% pbase et b=O => (t(pbase) %*% pbase) ^(-1)
  # matrice pbase n'est pas carree, le calcul permet d'avoir l'inverse dans projecteur
  # mais a priori comme orthonormee devrait être egale a t(pbase)
  # +/- le polynome constant
  # donc:
  projecteur=t(pbase)

  # projecteur de taille (degree+1) x ncol(simuls)

  ## calcul de H (H=Y.L^(-1)):
  # t(L) %*% t(projecteur) = I et projecteur %*% L = I
  # H*t(L) = Y => H*t(L) * t(projecteur) = Y * t(projecteur) => H =  Y * t(projecteur)

  coefpoly = as.matrix(simuls)%*% t(projecteur);
  # coefpoly de taille nrow(simuls) x (degree+1)

  return(list(H=coefpoly, # taille nrow(simuls) x (degree+1)
              L=pbase,    # taille ncol(simuls) x (degree+1)
              call.info=list(reduction="polynomial")))

}

