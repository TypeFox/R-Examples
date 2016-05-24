create.CCM <-
function(X.test, X.train = NULL, method = "pearson", use = "everything", verbose = 1) {
  if (is.null(X.train)) {
	K = cor(X.test, method = method, use = use)
	diag(K) = NA
  	attr(K, "class") <- "CCM"
        return(K)
  }

  if (is.null(rownames(X.train)) | is.null(rownames(X.test))) {
	if (nrow(X.train) != nrow(X.test)) {
	  cat("error: X.test and X.train must have the same number of rows or matching rownames must exist\n")
	  return(NULL)
	}
	m = 1:nrow(X.train)
  }
  m = match(rownames(X.train), rownames(X.test))
  if (method == "spearman") {
        if (verbose) cat("calculating ranks...\n")
  	X.test = apply(X.test,2,rank)
        X.train = apply(X.train,2,rank)
  }
  if (verbose) cat("calculating cors...\n")
  K = matrix(NA, nrow = ncol(X.test[m,]), ncol = ncol(X.train))
  for (i in 1:ncol(X.test)) {
	K[i,] = apply(X.train,2,cor,X.test[m,i])
  }
  attr(K, "class") <- "CCM"
  return(K)
}

