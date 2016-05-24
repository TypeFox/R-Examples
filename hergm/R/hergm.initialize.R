###########################################################################
# Copyright 2009 Nobody                                                   #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

hergm.initialize <- function(network, k, perturb)
# input: network, number of blocks
# output: membership indicators obtained by spectral clustering
{
  # Extract number of nodes:
  n <- network$gal$n

  # Extract adjacency matrix from network
  m <- as.matrix.network(network, matrix.type="adjacency")
  
  # Check whether network is symmetric:
  if (!isSymmetric(m)) 
    {
    # If not symmetric, make network symmetric by using OR rule:
    for (i in 1:n) 
      {
      for (j in 1:i) 
        {
        if (network[i,j]==1 || network[j,i]==1) 
          {
	  network[i,j] = network[j,i] = 1
	  }
        }
      } 
    } 

  # Spectral decomposition of network:
  b <- eigen(m, symmetric=TRUE)
  
  # k-means clustering of nodes based on "spectral.decomposition":
  c <- kmeans(b$vectors, centers=k, nstart=10) 

  # Extract membership indicators from "clustering":
  c.indicators <- c$cluster

  # Write membership indicators to console:
  cat("Membership Indicators:", c.indicators)
  
  if (perturb) 
    {
    # Sample membership indicators from S-neighborhood of vector of membership indicators:
    # - sampling S from uniform distribution on {1, ..., n/3}, where n is the number of nodes
    # - sampling S nodes without replacement
    # - move each of the S sampled nodes to another block by sampling from the uniform distribution on the other block indicators
    c.indicators <- sample(c.indicators, replace=TRUE)
    
    # Write sampled membership indicators to console:
    cat("\n")
    cat("Sampled Membership Indicators:", c.indicators)
    cat("\n")
    }

  # Return sampled partition, labeling blocks 0..k-1:
  indicator <- c.indicators - 1 # Replace by sampled membership indicators

}

