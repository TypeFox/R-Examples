# Multisensi R package ; file basis.mine.r (last modified: 2016-02-03) 
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
basis.mine <- function(simuls, basis.args=list(baseL=1*outer(sort(0:(ncol(simuls)-1)%%5),0:4,"==")))
#===========================================================================
{
  ## fait la decomposition sur une base donnee dans les arguments

  ## INPUTS
  ## simuls     : data.frame sur laquelle faire la decomposition (sortie des simulations du modele)
  ## basis.args : arguments specifiques a la methode de decomposition sous forme de list
  ## - baseL    : matrice des vecteurs de la base de decomposition
  ##              doit etre de taille ncol(simuls) x nb_vecteurs_de_la_base

  ## OUTPUTS
  ## H          : coefficients associes a la base de decomposition
  ##              de taille nrow(simuls) x ncol(baseL) 
  ## L          : matrice des vecteurs de la base, cad argument baseL
  ##              de taille ncol(simuls) x ncol(baseL)
  ##  NB : simuls ~ H %*% t(L)
  ## call.info  : contient des indications sur la fonction utilisee, call.info$reduction="matrix"

  b.args=basis.args
  # baseL doit etre de taille ncol(simuls) x (dim)

  # definition des valeurs par defaut des arguments
  if(is.null(basis.args$baseL)){
    stop("Warning : basis.mine argument 'baseL' not given,\n")
  }else if(!is.matrix(basis.args$baseL)){stop("baseL is not a matrix")}
  if(nrow(basis.args$baseL)!=ncol(simuls)){
    stop("Warning : matrix baseL dimension not compatible with outputs, nrow(baseL) must be ",ncol(simuls),"....\n")
  }

  # test orthogonalite -> warning
  if(all(t(basis.args$baseL)%*%basis.args$baseL-diag(1,ncol(basis.args$baseL),ncol(basis.args$baseL))<=1e-10)){
    # t(baseL)%*%baseL == I
    # a priori dans la cas base orthonormee inverse devrait être egale a t(pbase)
    projecteur=t(basis.args$baseL)
  }else{
    cat("Be careful : basis.mine argument 'baseL'is not orthogonal (t(baseL)%*%baseL != Id) \n")
    projecteur = solve( t(basis.args$baseL) %*% basis.args$baseL) %*% t(basis.args$baseL)
    # "solve" solves the equation ‘a %*% x = b’ for ‘x’ 
    # ici a =t(base) %*% base et b=O => (t(base) %*% base) ^(-1)
    # matrice base n'est pas carree, le calcul permet d'avoir l'inverse dans projecteur
    # test norme -> warning 
    #if(!all(sqrt(colSums(basis.args$baseL^2)) ==1)){cat("Warning : basis.mine argument 'baseL'pas normee \n")} 
  }

  # projecteur de taille (dim) x ncol(simuls)

  ## calcul de H (H=Y.L^(-1)):
  # t(L) %*% t(projecteur) = I et projecteur %*% L = I
  # H*t(L) = Y => H*t(L) * t(projecteur) = Y * t(projecteur) => H =  Y * t(projecteur)

  coefH = as.matrix(simuls)%*% t(projecteur);
  # coefH de taille nrow(simuls) x (dim)

  return(list(H=coefH,             # taille nrow(simuls) x (dim)
              L=basis.args$baseL, # taille ncol(simuls) x (dim)
              call.info=list(reduction="matrix")))

}

