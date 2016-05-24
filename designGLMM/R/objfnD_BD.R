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
objfnD_BD <- function(des,ntmt,blksz,sigb,sige,means,probs){
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
  cm2<-cm
  cm<-rbind(c(1,rep(0,ntmt)),cm)
  cm<-cm/rowSums(cm^2)
  # Calculate determinant information matrix for contrasts
  # Return negative determinant
  im<-cm %*% t(X) %*% oP %*% X %*% t(cm)
  # if(det(im)<1e-10) ret<-0 else ret<--det(im2)
  if(det(im)<1e-10) ret<-0 else {
    vc<-solve(im)[-1,-1]
    im2<-solve(vc)
    ret<--det(im2)
  }
  ret
}

