### z.transform.R  (2004-01-15)
###
###    Variance-Stabilizing Transformations of the Correlation Coefficient
###    
###
### Copyright 2003-04 Korbinian Strimmer
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


# Fisher's z-transform 
z.transform = function(r)
{
  return( atanh(r) )
}


# Hotelling's second-order transform
hotelling.transform = function(r, kappa)
{
  z = z.transform(r)
  
  return( z - (3*z+r)/(4*kappa) )
}
