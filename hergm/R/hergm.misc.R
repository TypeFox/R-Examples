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

length_mcmc <- function(d1, d2, k, n, predictions)
# input: number of ergm terms d1, number of hergm terms d2, number of blocks k, number of nodes n
# output: number of elements stored on iteration of MCMC algorithm
{
  terms <- d1 + (2 * d2) + (d2 * (k + 1)) 
  if (d2 > 0) terms <- terms + n + k + k + 1 
  if (predictions == TRUE) terms <- terms + d1 + d2
  #terms <- d1 # Number of ergm terms
  #       + (2 * d2) # Number of mean and precision parameters of Gaussian baseline distribution of Dirichlet / stick-breaking prior
  #       + ((d2 + 1) * k) # Number of hergm terms
  #       + n # Number of category indicators  
  #       + k # Number of category sizes 
  #       + k # Number of category probabilities
  #       + 1 # Clustering parameter
  #       + d1 + d2 # Number of posterior predictions
  #print("terms")
  #print(terms)
  terms
}

number_clusters <- function(alpha, n)
# input: scaling parameter of Dirichlet process <alpha>, number of units <n>
# output: first two moments of number of non-empty clusters
# see: Teh (2010), Encyclopedia of Machine Learning
{
  mean <- alpha * (digamma(alpha + n) - digamma(alpha))
  variance <- mean + (alpha * alpha) * (trigamma(alpha + n) - trigamma(alpha))
  moments <- list()
  moments$mean <- mean
  moments$variance <- variance
  moments
}

bernoulli_map_mean_to_natural <- function(n, mu, directed)
# input: number of nodes <n>, mean-value parameter of Bernoulli random graph model <mu>, indicator of whether network is directed
# output: natural parameter of Bernoulli random graph model
{
  if (directed == FALSE) df <- n * (n - 1) / 2
  else df <- n * (n - 1)
  theta <- - log((df - mu) / mu)
  theta
}

