# Multisensi R package ; file fctspline.jbrf.r (last modified: 2015-10-19) 
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
## fonctions Bsplines de l'annexe 14 du mémoire de stage de J.Baudet (2009_09)

####################################################################
#===========================================================================
bspline <- function(x=seq(0,1,len=101),k = knots, i = 1 ,m=2)
#===========================================================================
{
# evaluate ith b-spline basis function of order m at the values in x, given knot locations in k
# by Simon Wood
  if ( m== -1 ) # base of recursion
  {
    res = as.numeric( x < k[i+1] & x >= k[i] )
  }
  else{
    if( k[i+m+1] - k[i] == 0) { res0 = 0 }
    else{
      z0 = ( x - k[i] ) / ( k[i+m+1] - k[i])
      res0 = z0 * bspline(x,k,i,m-1)
    }
    if( k[i+m+2] - k[i+1] == 0 ){ res1 = 0 }
    else{ 
      z1 = ( k[i+m+2] - x ) / ( k[i+m+2] - k[i+1] )
      res1 = z1 * bspline(x,k,i+1,m-1)
    }
    res = res0 + res1
  }
  res
}

####################################################################
# nouvelles fonctions
####################################################################
#===========================================================================
sesBsplinesNORM <-function(x=seq(0,1,len=101),knots = 5, m=2)
#===========================================================================
{
  if( length(knots) == 1) {
    k = c(rep(min(x),m+1), seq(min(x),max(x),length=knots), rep(max(x),m+1))
    nbases = knots + m }
  else {
    k = c(rep(min(knots),m+1), sort(unique(knots)), rep(max(knots),m+1))
    nbases = length(unique(knots)) + m }
  courbes = matrix( NA, nrow=length(x), ncol=nbases )
  for (jj in 1 : nbases) courbes[,jj] = bspline(x,k,jj,m)
  courbes[length(x) , nbases ] = 1
  #TEST : normons les courbes :
  # pour que apply(GSI.spl$L^2,2,sum)==1 cad colSums(GSI.spl$L^2) ==1
  normcourbes=sqrt(colSums(courbes^2)) 
  for (i in 1:nbases){
    courbes[,i]=courbes[,i]/normcourbes[i]
  }
  projecteur = solve( t(courbes) %*% courbes) %*% t(courbes)
  # "solve" solves the equation ‘a %*% x = b’ for ‘x’ 
  # ici a =t(courbes) %*% courbes et b=O => (t(courbes) %*% courbes) ^(-1)
  # matrice courbes n'est pas carree, le calcul permet d'avoir l'inverse dans projecteur
    colnames(courbes)=paste("B",1:ncol(courbes),sep="")
  list(x = x, bsplines = courbes, knots = k, projecteur = projecteur)
}
####################################################################
#===========================================================================
sesBsplinesORTHONORM <-function(x=seq(0,1,len=101),knots = 5, m=2)
#===========================================================================
{
  if( length(knots) == 1) {
    k = c(rep(min(x),m+1), seq(min(x),max(x),length=knots), rep(max(x),m+1))
    nbases = knots + m }
  else {
    k = c(rep(min(knots),m+1), sort(unique(knots)), rep(max(knots),m+1))
    nbases = length(unique(knots)) + m }
  courbes = matrix( NA, nrow=length(x), ncol=nbases )
  for (jj in 1 : nbases) courbes[,jj] = bspline(x,k,jj,m)
  courbes[length(x) , nbases ] = 1
  #Normons les courbes :
  # pour que apply(GSI.spl$L^2,2,sum)==1 cad colSums(GSI.spl$L^2) ==1
  normcourbes=sqrt(colSums(courbes^2)) 
  for (i in 1:nbases){
    courbes[,i]=courbes[,i]/normcourbes[i]
  }
  # on applique transformation pour que la base soit orthogonale
  # fonction qr de R
  tmp=qr(courbes)
  ocourbes=qr.Q(tmp) # t(ocourbes)%*%ocourbes == I

  projecteur = t(ocourbes) #solve( t(ocourbes) %*% ocourbes) %*% t(ocourbes)
  # "solve" solves the equation ‘a %*% x = b’ for ‘x’ 
  # ici a =t(ocourbes) %*% ocourbes et b=O => (t(ocourbes) %*% ocourbes) ^(-1)
  # matrice ocourbes n'est pas carree, le calcul permet d'avoir l'inverse dans projecteur
  colnames(ocourbes)=paste("O",1:ncol(courbes),sep="")
  list(x = x, osplines = ocourbes, knots = k, projecteur = projecteur)
}

