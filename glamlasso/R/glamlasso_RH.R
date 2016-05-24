#
#     Description of this R script:
#     Rotated H-transform of an array A by a matrix M.
#
#     Intended for use with R.
#     Copyright (C) 2015 Adam Lund
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

# H-transform of an array A by a matrix X 
H<-function(M, A){
  d <- dim(A)
  Amat <- matrix(A, nrow = d[1])
  MAmat <- M %*% Amat 
  array(MAmat, c(nrow(MAmat), d[-1]))
}

# Rotation of an array A 
Rotate <- function(A){
  d <- 1:length(dim(A)) 
  d1 <- c(d[-1], d[1]) 
  aperm(A, d1)
}

#' @name RH
#' 
#' @aliases glamlasso_RH Rotate H
#' 
#' @title The Rotated H-transform of a 3d Array by a Matrix
#' 
#' @description  This function is an implementation of the \eqn{\rho}-operator found in 
#' \cite{Currie et al 2006}. It forms the basis of the GLAM arithmetic. 
#' 
#' @details For details see \cite{Currie et al 2006}. Note that this particular implementation 
#' is not used in the optimization routines underlying the glamlasso procedure.
#' 
#' @usage RH(M, A)
#' 
#' @param M a \eqn{n \times p_1} matrix.
#' @param A a 3d array of size \eqn{p_1 \times p_2 \times p_3}. 
#' 
#' @return A 3d array of size \eqn{p_2 \times p_3 \times n}.
#' 
#' @references 
#' Currie, I. D., M. Durban, and P. H. C. Eilers (2006). Generalized linear
#' array models with applications to multidimensional
#' smoothing. \emph{Journal of the Royal Statistical Society. Series B}. 68, 259-280.
#' 
#' @author Adam Lund

RH <- function(M, A){
  
  Rotate(H(M, A))

}
