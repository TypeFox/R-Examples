## Copyright (C) 1999 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## usage: w = gausswin(n, a)
##
## Generate an n-point gaussian window of the given width. Use larger a
## for a narrow window.  Use larger n for a smoother curve. 
##
##     w = exp( -(a*x)^2/2 )
##
## for x = linspace(-(n-1)/n, (n-1)/n, n)

gausswin <- function( n, w = 2.5)  {
  exp( -0.5 * ( w/n * seq(-(n-1), n-1, by = 2) ) ^ 2 )
} 
