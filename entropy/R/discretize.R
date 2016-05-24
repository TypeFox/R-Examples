### discretize.R  (2013-07-15)
###
###    Discretize a continuous random variable
###
### Copyright 2013 Korbinian Strimmer
###
###
### This file is part of the `entropy' library for R and related languages.
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


discretize = function( x, numBins, r=range(x) )
{
  b = seq(from=r[1], to=r[2], length.out=numBins+1 )
  y = table( cut(x, breaks=b , include.lowest=TRUE) )

  return( y )
}

discretize2d = function( x1, x2, numBins1, numBins2, r1=range(x1), r2=range(x2) )
{
  b1 = seq(from=r1[1], to=r1[2], length.out=numBins1+1 )
  b2 = seq(from=r2[1], to=r2[2], length.out=numBins2+1 )

  y2d = table( cut(x1, breaks=b1, include.lowest=TRUE ), 
               cut(x2, breaks=b2, include.lowest=TRUE ) )

  return( y2d )
}



