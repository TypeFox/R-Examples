checkNames <- function(X1, X2, X1.name="X1", X2.name="X2") {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  d1 <- ncol(X1)
  d2 <- ncol(X2)
  
  if (d1!=d2) stop(X1.name, " and ", X2.name, " must have the same numbers of columns")
  d <- d1
  
	
  # check names
  nm1 <- unique(colnames(X1))
  if (sum(nchar(nm1))==0) stop(X1.name, " does not contain any names")
  if (length(nm1)<d) stop("not enough names (ties?) found in ", X1.name)
  
  nm2 <- unique(colnames(X2))
  if (sum(nchar(nm2))==0) {
    warning(X2.name, " not named: ", X2.name, "'s variables are inherited from ", X1.name)
    colnames(X2) <- colnames(X1)
  } else {
    if (length(nm2)<d) stop("not enough names (ties?) found in ", X2.name)
    if (!all(nm2 %in% nm1)) stop("one name in ", X2.name, " is not in ", X1.name)
    if (!all(nm1 %in% nm2)) stop("one name in ", X1.name, " is not in ", X2.name)
    ind <- 1:length(nm1)
    names(ind) <- nm2
    X2 <- X2[, ind[nm1], drop=FALSE]
  }
  return(X2)
}