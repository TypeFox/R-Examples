### mpower.R  (2010-01-15)
###
###    Compute the Power of a Real Symmetric Matrix
###
### Copyright 2008-10 Korbinian Strimmer
###
###
### This file is part of the `corpcor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


# compute m^alpha where m is a symmetric matrix
mpower = function(m, alpha, pseudo=FALSE, tol)
{
    if( any( abs(m-t(m)) > 100*.Machine$double.eps  ) ) 
      stop("Input matrix is not symmetric!")

    em = eigen(m, symmetric = TRUE)
    eval = em$values

    # set small eigenvalues to exactly zero
    if( missing(tol) )
        tol = max(dim(m))*max(abs(eval))*.Machine$double.eps
    eval[abs(eval) <= tol] = 0

    if (pseudo) # use only the nonzero eigenvalues
    {
        idx = (eval != 0)
    }
    else # use all eigenvalues
    {
        idx = (1:length(eval))
    }

    e2 = eval[idx]^alpha
    ma = em$vectors[,idx, drop=FALSE] %*% 
         tcrossprod(diag(e2, nrow=length(e2)), 
         em$vectors[,idx, drop=FALSE])

    rownames(ma) = rownames(m)
    colnames(ma) = colnames(m)

    return(ma)
}


