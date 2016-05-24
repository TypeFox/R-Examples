
"null.space" <- function(x)
  { tmp <- svd(x)
    return(tmp$v[,which(tmp$d!=0)])
  }
