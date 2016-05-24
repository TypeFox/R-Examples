### kappa2n.R (2004-01-20)
###
###    Conversion of kappa to n and vice versaa 
###    
###
### Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
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



# sample size corresponding to kappa and G
kappa2n = function(kappa, p=2)
{
  return( kappa+p-1 )
}

# sample size corresponding to N and G
n2kappa = function(n, p=2)
{
  return( n-p+1 )
}


