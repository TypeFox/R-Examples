compute_norms_cat = function(X, res, n, p, numLevels, numCores=1){
  .Call("R_compute_norms_cat", X, res, n, p, numLevels, numCores, double(p))
}

compute_norms_cont = function(Z, res, n){
  res = matrix(res, nrow=1)
  abs(res %*% Z) / n
}

compute_norms_cat_cat = function(X, res, n, numLevels, indices, numCores=1){
  numInt = nrow(indices)
  .Call("R_compute_norms_cat_cat", X, res, n, numInt, numLevels, indices[, 1], indices[, 2], numCores, double(numInt))
}

compute_norms_cat_cont = function(X, Z, catNorms, res, n, numLevels, indices, numCores=1){
  numInt = nrow(indices)
  .Call("R_compute_norms_cat_cont", X, Z, catNorms, res, n, numInt, numLevels, indices[, 1], indices[, 2], numCores, double(numInt))
}

compute_norms_cont_cont = function(Z, contNorms, res, n, indices, numCores=1){
  numInt = nrow(indices)
  .Call("R_compute_norms_cont_cont", Z, contNorms, res, n, numInt, indices[, 1], indices[, 2], numCores, double(numInt))
}

