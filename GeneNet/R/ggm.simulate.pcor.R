### ggm.simulate.pcor  (2004-09-15)
###
###     Simulate GGM Networks
###
### Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
###
###
### This file is part of the `GeneNet' library for R and related languages.
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


# chose random GGM network
#  returns a positive definite partial correlation matrix
#  with the given proportion of non-zero elements
ggm.simulate.pcor = function(num.nodes, etaA=0.05)
{
  # num.nodes:  number of nodes in the network
  # etaA:          proportion of edges in the graph
  
  ####################
  eps = 0.0001
  # eps:        additional additive component to diag(precision_matrix)
  ####################
  
  
  num.edges = num.nodes*(num.nodes-1)/2
  num.elements = ceiling(num.edges*etaA) 

  # determine position of the non-zero elements
  element.idx = sample(1:num.edges, num.elements)
  precision.lo = rep(0,num.edges)
  
  # draw number 
  precision.lo[element.idx] = runif(num.elements,-1.0,+1.0)
 
  # construct symmetric matrix
  precision = matrix(0, nrow=num.nodes, ncol=num.nodes)
  precision[lower.tri(precision)] = precision.lo
  for(i in 1:(num.nodes-1))
  {
    for(j in (i+1):num.nodes)
    {
	precision[i,j] = precision[j,i]
    }
  }
	
  # construct diagonally dominant matrix (so that it is positive definite)	
  for(i in 1:num.nodes)
  {
    diag(precision)[i] = sum(abs(precision[,i]))+eps	
  }	
  pcor = cov2cor(precision)

 
  return(pcor)
}
