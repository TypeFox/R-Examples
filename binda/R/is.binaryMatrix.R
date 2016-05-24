### is.binaryMatrix.R  (2015-02-28)
###
###    Test For Binary Matrix
###
### Copyright 2013-15  Sebastian Gibb and Korbinian Strimmer
###
###
### This file is part of the `binda' library for R and related languages.
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

is.binaryMatrix = function(m)
{
  # is.matrix(m) && all(m == 1L | m == 0L) # slow
  is.matrix(m) && all(m*(m - 1) == 0L)
}

