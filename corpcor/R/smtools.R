### smtools.R  (2006-06-02)
###
###     Convert symmetric matrix to vector and back
###
### Copyright 2003-06 Korbinian Strimmer
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


# convert symmetric matrix to vector
sm2vec = function(m, diag = FALSE)
{
    return( as.vector(m[lower.tri(m, diag)]) )
}

# corresponding indices
sm.index = function(m, diag = FALSE)
{
  m.dim = length(diag(m))
 
  if (diag == TRUE)
    num.entries = m.dim*(m.dim+1)/2
  else
    num.entries = m.dim*(m.dim-1)/2
    
  index1 = rep(NA, num.entries )
  index2 = rep(NA, num.entries )

  if (diag == TRUE)
    delta = 0
  else
    delta = 1

  z = 1
  for (i in 1:(m.dim-delta))
    for (j in (i+delta):m.dim)
    {
      index1[z] = i
      index2[z] = j
      z = z+1
    }
      
 return( cbind(index1, index2) )
}

# convert vector to symmetric matrix
#
# note: if diag=FALSE then the diagonal will consist of NAs
#
vec2sm = function(vec, diag = FALSE, order = NULL)
{
  # dimension of matrix
  n = (sqrt(1+8*length(vec))+1)/2
  if (diag == TRUE) n = n-1
  if ( ceiling(n) != floor(n) )
    stop("Length of vector incompatible with symmetric matrix")
       
  # fill lower triangle of matrix     
  m = matrix(NA, nrow=n, ncol=n)
  lo = lower.tri(m, diag)
  if (is.null(order))
  {
    m[lo] = vec
  }
  else
  {
    # sort vector according to order
    vec.in.order = rep(NA, length(order))
    vec.in.order[order] = vec
    m[lo] = vec.in.order
  }
  
  # symmetrize
  for (i in 1:(n-1))
    for (j in (i+1):n)
         m[i, j] = m[j, i]   
  
  return( m )
}
