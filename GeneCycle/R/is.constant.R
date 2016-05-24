### is.constant.R  (2006-02-15)
###
###     Simple check for constant time series
###
### Copyright 2003-06 Korbinian Strimmer
###
###
### This file is part of the `GeneCycle' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
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



# checks wether a vector (time series) is constant
is.constant.single = function (v)
{
    return( all(diff(v) == 0))
}


# dito, but also for an array of time series
is.constant <- function(x)
{
  return( apply(as.matrix(x), 2, is.constant.single) )    
} 
