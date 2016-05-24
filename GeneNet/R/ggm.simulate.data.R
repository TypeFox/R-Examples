### ggm.simulate.data  (2013-06-16)
###
###     Simulate GGM data
###
### Copyright 2003-13 Juliane Schaefer and Korbinian Strimmer
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


# Simulate data from a given GGM model
#  input:  matrix with partial correlations 
#  output: multinormal data (with mu=0 and var=1)
ggm.simulate.data = function(sample.size, pcor)
{
    mu = rep(0, dim(pcor)[1])
  
    cor.mat = pcor2cor(pcor)
  
    return( myrmvnorm(sample.size, mu, cor.mat) )
}


########## internal ##########

# modified from mvtnorm package
# generate multinormal data with given mean vector and covariance 
myrmvnorm = function(n, mean, sigma)
{
  ev = eigen(sigma, symmetric = TRUE)
  tmp = ev$vectors %*% ( t(ev$vectors) * sqrt(ev$values) )
  tmp = matrix(rnorm(n * ncol(sigma)), nrow = n) %*% tmp
  tmp = sweep(tmp, 2, mean, "+")

  return(tmp)
}

