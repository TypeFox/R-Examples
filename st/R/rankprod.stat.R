### rankprod.R  (2008-12-18)
###
###    Two-sided rank products statistic
###
### Copyright 2008 Korbinian Strimmer
###
### This function is in part based on code from Henry Wirth.
###
### This file is part of the `st' library for R and related languages.
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


rankprod.stat = function (X, L)
{
  FUN = rankprod.fun(L=L)
  score = FUN(X)
  
  return( score )
}

rankprod.fun <- function (L)
{
  if (missing(L)) stop("Class labels are missing!")
  L = factor(L)
  cl = levels(L)
  if (length(cl) != 2) stop("Class labels must be specified for two groups, not more or less!")

  idx1 = which( L == levels(L)[1] )
  idx2 = which( L == levels(L)[2] )

  function(X)
  {
    ranks = array(0, dim=c( ncol(X), length(idx1)*length(idx2) ) )

    for( i in 1:length(idx1) )
    {
      for( j in 1:length(idx2) )
      {
        index = j+length(idx2)*(i-1)   
        ranks[ , index] = rank( -abs( X[idx1[i], ] - X[idx2[j], ] ) )
      }
    } 
       
    # geometric mean         
    avg.rank = exp(rowMeans(log(ranks)))

    # return rank complement
    compl.rank = ncol(X)-avg.rank+1

    return( compl.rank ) 
  }
}

