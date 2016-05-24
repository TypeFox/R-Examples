# Scale importance measures
# imp: matrix of predictors x importance type
# impSD: vector of standard deviation of importances for predictors
.scaleImp <- function(imp, impSD) {
  imp[, -ncol(imp)] <- imp[, -ncol(imp), drop = FALSE] / 
    ifelse(impSD < .Machine$double.eps, 1, impSD)
  return(imp)
}

# Make a 3-dimensional array of importance values
# x: initial 2-dimensional array of predictors x importance type
# z.dim: dimensionality of z axis
# z.dimnames: dimnames of z axis
.makeImpArray <- function(x, z.dim, z.dimnames) {
  array(
    x, dim = c(nrow(x), ncol(x), z.dim), 
    dimnames = list(rownames(x), colnames(x), z.dimnames)
  )
}

# Calculate importance p-values from null distribution and observed importances
# rp: rfPermute object
#' @importFrom swfscMisc pVal
.calcImpPval <- function(rp) {
  calcPval <- function(obs, null) {
    t(sapply(1:nrow(null), function(i) {
      sapply(1:ncol(null), function(j) pVal(obs[i, j], null[i, j, ]))
    }))
  }
  imp <- rp$importance
  impSD <- rp$importanceSD
  null.dist <- rp$null.dist
  rm(rp)
  arr <- .makeImpArray(calcPval(imp, null.dist$unscaled), 2, NULL)
  arr[, , 2] <- calcPval(.scaleImp(imp, impSD), null.dist$scaled)
  dimnames(arr) <- list(rownames(imp), colnames(imp), c("unscaled", "scaled"))
  return(arr)
}