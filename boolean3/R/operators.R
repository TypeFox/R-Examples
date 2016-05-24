### ----------------------------------------------------------------------------
### This file is part of boolean3

### Copyright (C) 2011--2014 Jason W. Morgan <morgan.746@osu.edu>

### boolean3 represents a substantial re-write of the original boolean package
### developed by Bear Braumoeller, Ben Goodrich, and Jacob Kline. This version
### was developed under the direction of Bear Braumoeller and with support from
### The Ohio State University's College of Social and Behavioral Sciences.

### boolean3 and is free software: you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the Free
### Software Foundation, either version 3 of the License, or (at your option)
### any later version.

### This program is distributed in the hope that it will be useful, but WITHOUT
### ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
### FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
### more details.

### You should have received a copy of the GNU General Public License along with
### this program.  If not, see <http://www.gnu.org/licenses/>.

### ----------------------------------------------------------------------------


## Boolean `or' operator.
## @name or
## @param x vector of type numeric.  
## @param y vector of type numeric.
## @return vector of type numeric.
## @author Jason Morgan
or  <- function(x, y) {
  ##1 - (1 - x) * (1 - y)
  ## This version below is slightly faster as it performs one less calculation.
  x + y - x*y
}

## Boolean `and' operator.
##
## Boolean `and' operator.
## @name and
## @param x vector of type numeric.  
## @param y vector of type numeric.
## @return vector of type numeric.
## @author Jason Morgan
and <- function(x, y) {
  x * y
}

"%|%" <- function(x, y) {
  or(x,y)
}

"%&%" <- function(x, y) {
  and(x,y)
} 

