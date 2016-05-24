### rebuild.cov.R (2006-05-25)
###
###    Rebuild and Decompose (Inverse) Covariance Matrix
###    
###
### Copyright 2003-06 Korbinian Strimmer
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




# rebuild covariance matrix from correlations (r) and variances (v)
rebuild.cov = function(r, v)
{
  if ( any( v < 0) ) stop("Negative variance encountered!")
  
  sd = sqrt(v)
  m = sweep(sweep(r, 1, sd, "*"), 2, sd, "*") 
    
  return(m)
}

# decompose covariance matrix into correlations and variances
decompose.cov = function(m)
{
  v = diag(m)
  r = cov2cor(m)

  return( list(r=r, v=v) )
}

# rebuild precision matrix from partial correlations (pr) and partial variances (pv)
rebuild.invcov = function(pr, pv)
{
  if ( any( pv < 0) ) stop("Negative partial variance encountered!")
    
  ipsd = sqrt(1/pv)
  m = -sweep(sweep(pr, 1, ipsd, "*"), 2, ipsd, "*") 
  diag(m) = -diag(m)
  
  return(m)
}

# decompose precision matrix into partial correlations and partial variances
decompose.invcov = function(m)
{
  pv = 1/diag(m)
  
  m = -m
  diag(m) = -diag(m)
  pr = cov2cor(m)

  return( list(pr=pr, pv=pv) )
}








