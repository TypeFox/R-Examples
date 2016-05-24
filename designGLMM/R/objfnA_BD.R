# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Define an objective function
#' @export
objfnA_BD <- function(des,ntmt,blksz,sigb,sige,means,probs=c(1)){
  nblk<-length(des)/blksz             # Calculate number of blocks
  xmeans<- means[des]                 # Create a vector of conditional means for each unit
  diagall<-sigb^2/(sige^2+1/xmeans)       # Create a vector of variance ratios
  ellvec<-matrix(diagall,nrow=blksz)  # Partition diagall into blocks
  # Calculate the diagonal blocks of projection matrix
  blk<-apply(ellvec,2,function(x) diag(x)/sigb^2 - x%*%t(x)*1/(sigb^2*(1+sqrt(t(x)) %*% sqrt(x)))[1,1])
  oP<-matrix(0,length(des),length(des))  # Arrange into block diagonal matrix
  for(q in 1:nblk){
    oP[((q-1)*blksz+1):((q-1)*blksz+blksz),((q-1)*blksz+1):((q-1)*blksz+blksz)]<-matrix(blk[,q],nrow=blksz)
  }
  # Calculate design matrix X
  Xmaster<-cbind(rep(1,ntmt),diag(1,ntmt))
  X<-t(sapply(des,function(x) Xmaster[x,]))
  # Set up contrasts
  cm<-matrix(0,nrow=ntmt,ncol=ntmt)
  cm[upper.tri(cm, diag = FALSE)]<--1
  cm<-cm+diag((ntmt-1):0)
  cm<-cbind(rep(0,ntmt),cm)
  cm<-matrix(cm[1:(ntmt-1),],nrow=(ntmt-1))
  cm<-rbind(c(1,rep(0,ntmt)),cm)
  cm<-cm/rowSums(cm^2)
  # Calculate inverse of information matrix for contrasts (if possible)
  # Return trace of the inverse of the information matrix if system is identifiable and a large number otherwise
  inversematrix<-try(solve(cm %*% t(X) %*% oP %*% X %*% t(cm)),silent=TRUE)
  if(class(inversematrix)=="try-error") ret<-1e160 else ret<-sum(diag(inversematrix))-inversematrix[1,1]
  ret
}

