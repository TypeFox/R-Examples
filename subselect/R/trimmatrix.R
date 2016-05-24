trim.matrix <- function(mat,tolval=10*.Machine$double.eps)
{
  p <- dim(mat)[2]
  matindices <- 1:p                  
  mat.eig <- eigen(mat,symmetric=TRUE)
  discard <- rep(FALSE,p)
  newmat <- mat
  newmatindices <- matindices        
  while(mat.eig$values[p]/mat.eig$values[1] < tolval)
  {
    int <- as.numeric(newmatindices[order(abs(mat.eig$vectors[,p]),decreasing=TRUE)[1]])   
    discard[int] <- TRUE
    newmat <- mat[!discard,!discard]
    newmatindices <- matindices[!discard]   
    p <- p-1
    mat.eig <- eigen(newmat,symmetric=TRUE)   
  }
  size <- dim(newmat)[2]
  output <- list(newmat,as.numeric(matindices[discard]),colnames(mat)[discard],size)   
  names(output) <- c("trimmedmat","numbers.discarded","names.discarded","size")
  output
}