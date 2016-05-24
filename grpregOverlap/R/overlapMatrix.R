
## function: overlap matrix c[i,j] = # of overlaps between group i and j.
##           Diagonals are group size. Provide interface for user to check 
##           overlapping structure.
# ------------------------------------------------------------------------------
overlapMatrix <- function(X, group) {
  inc.mat <- incidenceMatrix(X, group)
  over.mat <- Matrix(inc.mat %*% t(inc.mat), sparse = TRUE, dimnames = dimnames(inc.mat))
  over.mat
}
# ------------------------------------------------------------------------------

## function: incidence matrix: I[i, j] = 1 if group i contains variable j.
# ------------------------------------------------------------------------------
incidenceMatrix <- function(X, group) {
  n <- nrow(X)
  p <- ncol(X)
  if (class(group) != 'list') {
    stop("Argument 'group' must be a list of integer indices or character names of variables!")
  }
  J <- length(group)
  grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE, 
                    dimnames=list(as.character(rep(NA, J)),
                                  as.character(rep(NA, p))))    
  if(is.null(colnames(X))) {
    colnames(X) <- paste("V", 1:ncol(X), sep="")    
  }
  if (is.null(names(group))) {
    names(group) <- paste("grp", 1:J, sep="")
  }
  
  if (class(group[[1]]) == 'numeric') {
    for (i in 1:J) {
      ind <- group[[i]]
      grp.mat[i, ind] <- 1
      colnames(grp.mat)[ind] <- colnames(X)[ind]
    }
  } else { ## character, names of variables
    for (i in 1:J) {
      grp.i <- as.character(group[[i]])
      ind <- colnames(X) %in% grp.i
      grp.mat[i, ] <- 1*ind
      colnames(grp.mat)[ind] <- colnames(X)[ind]
    }
  }
  rownames(grp.mat) <- as.character(names(group))
  # check grp.mat
  if (all(grp.mat == 0)) {
    stop("The names of variables in X don't match with names in group!")
  }
  
  ## TODO:
  ## (1) handle cases where variables not belong to any of groups in 'group'
  ## put each variable into a separate group, stack those group at right
  ## (2) provide option of removing groups including only one variable.
  ## Will add this later...
  
#   grp.mat <- grp.mat[, colSums(grp.mat) != 0]
#   if (ncol(grp.mat) < p) { 
#     # there exist variables that don't belong to 'group'. 
#     # create a bottom-right corner identity matrix.
#     colnames.new <- c(colnames(grp.mat), setdiff(colnames(X), colnames(grp.mat)))
#     rownames.new <- c(rownames(grp.mat), 
#                       paste("grp", (J+1):(J+p-ncol(grp.mat)), sep=""))
#     grp.mat <- bdiag(grp.mat, Diagonal(p-ncol(grp.mat)))
#     dimnames(grp.mat) <- list(rownames.new, colnames.new)
#   }

  grp.mat
}
# ------------------------------------------------------------------------------
