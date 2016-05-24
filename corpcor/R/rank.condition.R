### condition.R  (2013-5-15)
###
###     Rank, condition, and positive definiteness of a matrix
###
### Copyright 2003-13 Korbinian Strimmer
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


# checks whether a matrix is positive definite
is.positive.definite = function (m, tol, method=c("eigen", "chol"))
{
    method = match.arg(method)
    
    if (!is.matrix(m)) m = as.matrix(m)

    if (method=="eigen")
    {
        eval = eigen(m, only.values = TRUE, symmetric=TRUE)$values
        if (is.complex( eval ))
        {
           warning("Input matrix has complex eigenvalues!")
           return(FALSE)
        }

        if( missing(tol) )
            tol = max(dim(m))*max(abs(eval))*.Machine$double.eps
   
        if (sum(eval > tol) == length(eval))
            return(TRUE)
        else
            return(FALSE)
    }
    
    if (method=="chol")
    {
	val = try(chol(m), silent=TRUE)
  
        if (class(val) == "try-error")
            return(FALSE)
        else
            return(TRUE)    
    }
}


# Method by Higham 1988
make.positive.definite = function(m, tol)
{
  if (!is.matrix(m)) m = as.matrix(m)

  d = dim(m)[1] 
  if ( dim(m)[2] != d ) stop("Input matrix is not square!")
   
  es = eigen(m, symmetric=TRUE)
  esv = es$values
  
  if (missing(tol))
      tol = d*max(abs(esv))*.Machine$double.eps 
  delta =  2*tol # factor to is just to make sure the resulting
                  # matrix passes all numerical tests of positive definiteness
  
  tau = pmax(0, delta - esv)
  dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)    
  
  #print(max(DA))
  #print(esv[1]/delta)
      
  return( m +  dm )
}




# rank and condition of a matrix 
rank.condition = function (m, tol)
{
    d = svd(m, nv=0, nu=0)$d # compute only singular values
    
    max.d = d[1]
    min.d = d[length(d)]
    
    if( missing(tol) ) 
        tol = max(dim(m))*max.d*.Machine$double.eps
    
    r = sum(d > tol) # rank: number of singular values larger than tol
    
    if (r < min(dim(m)) ) min.d = 0 # if matrix is singular then set the  smallest
                                     # singular value to 0, and hence condition = INF
    
    c = max.d/min.d
    
    return(list(rank = r, condition = c, tol=tol))
}

