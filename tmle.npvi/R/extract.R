extract <- function(mat, exclude) {
  theVar <- setdiff(colnames(mat), exclude)
  mat[, theVar, drop=FALSE]
}

extractW <- function(#Extracts W Columns from Matrix of Observations
### Extracts the \code{W} column(s) from a \code{matrix} of observations.
    mat
### A   \code{matrix}  of   observations,  as   the  \code{obs}   argument  of
### \code{function} \code{\link{tmle.npvi}}.
    ) {
  ##details<< Mainly for internal use.
  extract(mat, exclude=c("X", "Y"))
### The \code{matrix}  extracted from \code{mat} by removing  the two \code{X}
### and \code{Y} columns. 
}

extractXW <- function(#Removes the Y Column from Matrix of Observations
### Removes the \code{Y} column from a \code{matrix} of observations.
    mat
### A   \code{matrix}  of   observations,  as   the  \code{obs}   argument  of
### \code{function} \code{\link{tmle.npvi}}.    
    ) {
  ##details<< Mainly for internal use.
  mat <- extract(mat, exclude="Y")

  ## enforcing the order to be X, W
  cn <- colnames(mat)
  wX <- which(cn=="X")
  cbind(mat[, wX, drop=FALSE], mat[, -wX, drop=FALSE])
### The  \code{matrix}  extracted from  \code{mat}  by  removing the  \code{Y}
### column in such a way that the first column is \code{X}.
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

