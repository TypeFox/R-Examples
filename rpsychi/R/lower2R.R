lower2R <-
function(x, varname=NULL){
  temp <- function(y){
    (y ^2 - y)/2 - length(x)
  }

  if(length(x) != 1){
    ncol <- round(uniroot(temp, interval=c(0, length(x)))$root)
  }else{
    ncol <- 2
  }
  
  mat <- matrix(0, ncol=ncol, nrow=ncol)
  mat[lower.tri(mat)] <- c(x)
  mat <- mat + t(mat)
  diag(mat) <- 1
  
  if(is.null(varname)){
    varname <- "z"
    colnames(mat) <- rownames(mat) <- paste(varname,1:ncol(mat), sep="")
  }else{
    colnames(mat) <- rownames(mat) <- varname
  }
  
  return(mat)
}
