# =================================================================
# This function provides a formula interface to the function
# knn() in package cluster, that implements k-nearest neighbours
# classifiers.
# On top of that function this adds the possibility of normalizing
# the data through parameter norm (see examples below)
# -----------------------------------------------------------------
# Mar 2010, Luis Torgo
# -----------------------------------------------------------------
kNN <- function(form,train,test,norm=T,norm.stats=NULL,...) {
  ##require(class,quietly=TRUE)
  tgtCol <- which(colnames(train)==as.character(form[[2]]))
  if (norm) {
    if (is.null(norm.stats)) tmp <- scale(train[,-tgtCol],center=T,scale=T)
    else tmp <- scale(train[,-tgtCol],center=norm.stats[[1]],scale=norm.stats[[2]])
    train[,-tgtCol] <- tmp
    ms <- attr(tmp,"scaled:center")
    ss <- attr(tmp,"scaled:scale")
    test[,-tgtCol] <- scale(test[,-tgtCol],center=ms,scale=ss)
  }
  knn(train[,-tgtCol],test[,-tgtCol],train[,tgtCol],...)
}
