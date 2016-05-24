################################################################################
################################################################################
## Direct and indirect preferences

# Copyright 2014 Julien Boelaert.
# 
# This file is part of revealedPrefs.
# 
# revealedPrefs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# revealedPrefs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with revealedPrefs.  If not, see <http://www.gnu.org/licenses/>.

## Function to compute all direct preferences
## afriat.par in (0,1), 1 for standard RP
## prefs[i, j]= 0 iff i not prefered to j (afriat.par * p_i q_i < p_i q_j)
## prefs[i, j]= 1 iff equality with prices i (afriat.par * p_i q_i == p_i q_j)
## prefs[i, j]= 2 iff i strictly prefered to j (afriat.par * p_i q_i > p_i q_j)
directPrefs <- function(x, p, afriat.par= 1) {
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  px <- p %*% t(x)
  prefs <- matrix(0, nrow= nrow(x), ncol= nrow(x))
  prefs[afriat.par * diag(px) == px] <- 1
  prefs[afriat.par * diag(px) > px] <- 2
  
  prefs
}

## Function to compute all indirect preferences
## prefs[i, j]= 0 iff i not indirectly prefered to j
## prefs[i, j]= 1 iff i indirectly prefered to j (only equalities)
## prefs[i, j]= 2 iff i indirectly strictly prefered to j 
##                    (with at least one strict preference)
indirectPrefs <- function(x, p, afriat.par= 1) {
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  .Call("IndirectPrefs", p %*% t(x), afriat.par, PACKAGE= "revealedPrefs")
}
