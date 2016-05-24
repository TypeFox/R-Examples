##' @export
formula.lvm <- function(x,char=FALSE,all=FALSE,...) {
  A <- index(x)$A
  res <- c()
  for (i in seq_len(ncol(A))) {
    if (all || !(colnames(A)[i]%in%c(index(x)$exogenous,parameter(x)) )) {
      f <- paste(colnames(A)[i],"~ 1")
      if (any(A[,i]!=0)) {
        f <- (paste(colnames(A)[i],"~",paste(colnames(A)[A[,i]!=0],collapse="+")))
      }
      if (!char)
        f <- formula(f)
      res <- c(res, list(f))
    }
  }
  return(res)
}


##' @export
formula.lvmfit <- formula.lvm
