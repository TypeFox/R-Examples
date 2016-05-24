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
objfnA_CRD <- function(des,ntmt,sige,means,probs=c(1)){
  # Calculate oP matrix
  oP<-diag(1/(1/means[des]+sige^2))
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
  inversematrix<-try(solve(cm %*% t(X) %*% oP %*% X %*% t(cm)),silent=TRUE)
  if(class(inversematrix)=="try-error") ret<-1e160 else ret<-sum(diag(inversematrix))-inversematrix[1,1]
  ret
}

