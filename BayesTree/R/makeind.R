makeind <- function(x,all=TRUE){
  if(!is.data.frame(x)) cat('\n\twarning: x in makeind is not a dataframe\n')
  dummies.all <- NULL
  for (i in 1:ncol(x)) {
    if (is.factor(x[,i])){
      n.levels <- length(levels(x[,i]))
      dummies <- class.ind(x[,i])
      dimnames(dummies)[[2]] <- paste(dimnames(x)[[2]][i],dimnames(dummies)[[2]],sep='.')
      if ((n.levels==2) | !all) dummies <- dummies[,-ncol(dummies),drop=F]
      dummies.all <- cbind(dummies.all,dummies)
    }
  }
  if(!is.null(dummies.all)) {
    x <- x[,unlist(lapply(x,is.numeric))] # drop any var which is not numeric
    x <- as.matrix(cbind(x,dummies.all))  # add on all the dummies
  }
  return(as.matrix(x)) # coming in, x is a data.frame, coming out, it is a matrix
}
