###################################################################
# sivipm R package
# Copyright INRA 2016
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################

allmono <- function(p, deg) {
  # FUNCTION
  # List of all the monomials in a polynomial with p variables
  # of degree = deg
  # RETURN
  # An array  nmono x p.  For each monomial i, 
  # indic[i,j]= degree of the variable j or 0
  # AUTHOR
  # J.P. Gauchi, INRA/MaIAGE
  # EXAMPLE
  # > allmono(4,2)
  #         [,1] [,2] [,3] [,4]
  #  [1,]    0    0    0    2
  #  [2,]    0    0    1    1
  #  [3,]    0    0    2    0
  #  [4,]    0    1    0    1
  #  [5,]    0    1    1    0
  #  [6,]    0    2    0    0
  #  [7,]    1    0    0    1
  #  [8,]    1    0    1    0
  #  [9,]    1    1    0    0
  # [10,]    2    0    0    0
  # -------------------------------------------

  if (p<=1) stop("p should be greater than 1")
  if (deg<1) stop("deg should be greater than 0")
  valn <- p+deg-1
  lmono <- t(combn(valn, (p-1)))
  nmono <- nrow(lmono)
  indic <- matrix(NA, nrow=nmono, ncol=p)
  indic[,1] <- abs(lmono[,1]-1)
  clmono <- ncol(lmono)
  indic[,2:clmono] <- abs(lmono[,2:clmono] - lmono[,1:(clmono-1)] - 1)
  indic[,p] <- abs(lmono[, clmono] - valn)
  return(indic)
} # fin allmono
# -----------------------------------------------------
convertindic <- function(indic, deg) {
  # FUNCTION
  # Convert an indicatrice matrix into a list
  # INPUT
  # A matrix nmono X p
  # indic[i,j]=degree of the variable j in the monomial i
  # RETURN
  # A  list of nmono components.
  # For each monomial i, 
  # tab[[i]] is a vector of length = deg, which contains
  # the numbers of the variables in the monomial, completed with zeros
  # -------------------------------------------
  
  nmono <- nrow(indic)
  
  convert1row <- function (indic, deg) {
      liste <- c()
    for (j in 1:length(indic)) {
      if (indic[j] >0) {
        liste <- c(liste, rep(j, indic[j]))
      }
    } # fin j
    return(liste)
  } # fin convert1row
  tab <- apply(indic, 1, convert1row, deg)
  
  # tab is a list, because all the returns of convert1row
  # has not the same length
  # Fill in with zeros so each component of the list is
  # a vector of length deg
  completezero <- function(X, deg) {
         if (length(X) <deg) {
           X[(length(X)+1):deg] <- 0
         }
         X
       } # fin completezero
  tab <- lapply(tab, completezero, deg)
  return(tab)
}

# ---------------------------------------------------
polynomecomplet <- function(p, deg) {
  # FUNCTION
  # Generate a list coding for the full polynomial with p variables
  # of degree = deg, including all the monomials of degree <= deg 
  # RETURN
  # A  list of nmono components.
  # For each monomial i, 
  # tab[[i]] is a vector of length = deg, which contains
  # the numbers of the variables in the monomial, completed with zeros
  # EXAMPLE
  # > polynomecomplet(4,2)
  # [[1]]
  # [1] 4 0
  # 
  # [[2]]
  # [1] 3 0
  # ...
  # [[13]]
  # [1] 1 2
  # 
  # [[14]]
  # [1] 1 1
  # -------------------------------------------

  # 1/ generate an array  nmono x p
  res <- matrix(NA, nrow=0, ncol=p)
  for (d in 1:deg) {
    A <-  allmono(p, d)
    nlig <- nrow(A)
# Inversion of the lines to put the monomes in right order
    A[1:nlig,] <- A[nlig:1,]
    res <- rbind(res, A)
    }

  # 2/ convert into a nmono x deg list
  return(convertindic(res, deg))
} # fin polynomecomplet

