zeroInd <-
function(Amat, r){
  if (sum(t(Amat)!=Amat)>0){
    stop("This method only works for symmetric matrix!")
  }
  p <- dim(Amat)[1]
  oneMat <- matrix(0, p, p)
  zeroMat <- matrix(0, p, p)
  
  one.pos <- which(Amat!=0, arr.ind = TRUE)
  zero.pos <- which(Amat==0, arr.ind = TRUE)
  
  zero.pos <- zero.pos[which(zero.pos[,1] > zero.pos[,2]) ,]
  sel.zero <- sample(seq(1, dim(zero.pos)[1]), r * dim(zero.pos)[1], replace = FALSE) 
  zeroMat[zero.pos[sel.zero, ]] <- 1
  zeroMat <- zeroMat + t(zeroMat)  
  zeroArr <- zero.pos[sel.zero, ]
  
  out <- list()
  out$zeroArr = zeroArr
  out$zeroMat = zeroMat
  
  if (dim(one.pos)[1] == 0){
    warning("The matrix is zero!")
    out$oneMat = matrix(0, p, p)
  } else 
  {
    one.pos <- one.pos[which(one.pos[,1] > one.pos[,2]) ,]
    if (is.null(dim(one.pos))){
      one.pos = matrix(one.pos, nrow = 1)
    }
     
    sel.one <- sample(seq(1, dim(one.pos)[1]), r * dim(one.pos)[1], replace = FALSE) 
    oneMat[one.pos[sel.one, ]] <- 1
    oneMat <- oneMat + t(oneMat)
    diag(oneMat) <- 0
    
    out$oneMat = oneMat  
  }
  
  return(out)  
}
