initsvd <- function(X) {
  n = NROW(X)
  p = NCOL(X)
  ifelse(n>=p, return(svd(X,nu=0,nv=1)$v)
             , return(svd(X,nu=1,nv=0)$u)
        )
}
