### dyn.cor.R  (2008-11-14)
###
###    Dynamical Correlation and Covariance
###
### Copyright 2005-2008 Rainer Opgen-Rhein and Korbinian Strimmer
###
###
###
### This file is part of the `GeneTS' library for R and related languages.
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




######  correlation related ######

# dynamical correlation
dyn.cor = function(x, lambda, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  c = cor.shrink(x, lambda=lambda, w=w, verbose=verbose)

  return( c )
}

# dynamical inverse correlation
dyn.invcor = function(x, lambda, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  c = invcor.shrink(x, lambda=lambda, w=w, verbose=verbose)

  return( c )
}

# dynamical partial correlation
dyn.pcor = function(x, lambda, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  c = pcor.shrink(x, lambda=lambda, w=w, verbose=verbose)

  return( c )
}



######  covariance related ######

# dynamical covariance
dyn.cov = function(x, lambda, lambda.var, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  c = cov.shrink(x, lambda=lambda, lambda.var=lambda.var,
         w=w, verbose=verbose)

  return( c )
}


# dynamical inverse covariance
dyn.invcov = function(x, lambda, lambda.var, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  c = invcov.shrink(x, lambda=lambda, lambda.var=lambda.var, 
         w=w, verbose=verbose)

  return( c )
}


# dynamical variance
dyn.var = function(x, lambda.var, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  sv = var.shrink(x, lambda.var=lambda.var, 
          w=w, verbose=verbose)

  return( sv )
}


# dynamical partial variance
dyn.pvar = function(x, lambda, lambda.var, verbose=TRUE)
{
  w = dyn.weights(x)
  
  # weighted covariance
  pv = pvar.shrink(x, lambda=lambda, lambda.var=lambda.var, 
          w=w, verbose=verbose)

  return( pv )
}


