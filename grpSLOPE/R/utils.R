###############################################################################
#
#    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#' Get a groupID object
#'
#' Mostly intended for internal use.
#'
#' @param group A vector describing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#'
#' @return An object of class groupID, which is a list, whose members are 
#'    vectors of indices corresponding to each group. The names of
#'    the list members are the corresponding group names.
#'
#' @examples
#' group  <- c("A", "A", 2, 9, "A", 9, 9, 2, "A")
#' group.id <- getGroupID(group)
#' group.id
#' # $A
#' # [1] 1 2 5 9
#' # 
#' # $`2`
#' # [1] 3 8
#' # 
#' # $`9`
#' # [1] 4 6 7
#' # 
#' # attr(,"class")
#' # [1] "groupID"
#'
#' @export
getGroupID <- function(group) {
  group.unique <- unique(group)
  n.group <- length(group.unique)
  group.id <- list()
  for (i in 1:n.group){
    id <- as.character(group.unique[i])
    group.id[[id]] <- which(group==group.unique[i])
  }
  class(group.id) <- "groupID"
  return(group.id)
}

# Orthogonalize each group of columns of a matrix A.
# For i = 1, ..., m let A_i = A[ , group_i] and compute
# A_i[ , P] = Q %*% R, where P is a permutation vector.
#
orthogonalizeGroups <- function(X, group.id) {
  n.group <- length(group.id)

  getGroupQR <- function(ids) {
    submat <- X[ , ids]

    if (length(ids) == 1) {
      return(list(Q=as.matrix(submat), R=1, P=1))
    } else {
      submat.qr <- qr(submat, LAPACK=TRUE)
      return(list(Q=qr.Q(submat.qr),
                  R=qr.R(submat.qr),
                  P=submat.qr$pivot))
    }
  }

  return(lapply(group.id, getGroupQR))
}
