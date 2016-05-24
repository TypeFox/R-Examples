
corrmsrs <- function(X, minprop=0.1, maxprop=0.9) {
  if (is.vector(X)) X <- as.matrix(X)
  if (!is.matrix(X) && !is.data.frame(X)) 
    stop("X is not a matrix nor data.frame")
  if (minprop < 0 || maxprop > 1 || minprop > maxprop) 
    stop("minprop/maxprop are not correct proportions")
  ilist <- vector(ncol(X), mode="list")
  for (i in 1:ncol(X))
    ilist[[i]] <- irank(as.numeric(X[,i]))
  a <- .Call("newcorr", ilist=ilist, as.double(c(minprop, maxprop)),
             PACKAGE="maxstat")
  corrm <- a[[1]]
  coldel <- a[[2]]
  rowdel <- a[[3]]
  corrm <- corrm[rowdel == 0, coldel == 0]
  corrm
}
