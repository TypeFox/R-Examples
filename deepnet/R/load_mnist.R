##' Load MNIST DataSet
##'
##' Load MNIST DataSet
##' @param dir dir of minst dataset
##' @return mnist dataset
##'   train$n   number of train samples
##'   train$x   pix of every train sample image
##'   train$y   label of every train sample image
##'   train$yy  one-of-c vector of label of train sample image
##'   test$n   number of test samples
##'   test$x   pix of every test sample image
##'   test$y   label of every test sample image
##'   test$yy  one-of-c vector of label of test sample image
##' @author Xiao Rong
##' @export
load.mnist <- function(dir) {
  load.image.file <- function(filename) {
    ret <- list()
    f <- file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n <- readBin(f,'integer',n=1,size=4,endian='big')
    nrow <- readBin(f,'integer',n=1,size=4,endian='big')
    ncol <- readBin(f,'integer',n=1,size=4,endian='big')
    x <- readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x <- matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load.label.file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  mnist <- list()
  mnist$train <- load.image.file(paste(dir,'/train-images-idx3-ubyte',sep=""))
  mnist$test <- load.image.file(paste(dir,'/t10k-images-idx3-ubyte',sep=""))
  
  mnist$train$y <- load.label.file(paste(dir,'/train-labels-idx1-ubyte',sep=""))
  n <- length(mnist$train$y)
  mnist$train$yy <- matrix(rep(0,n*10),nrow=n,ncol=10)
  for (i in 1:n){
    mnist$train$yy[i,mnist$train$y[i] + 1] <- 1
  }
  mnist$test$y <- load.label.file(paste(dir,'/t10k-labels-idx1-ubyte',sep=""))
  m <- length(mnist$test$y)
  mnist$test$yy <- matrix(rep(0,m*10),nrow=m,ncol=10)
  for (j in 1:m){
    mnist$test$yy[j,mnist$test$y[j] + 1] <- 1
  }
  mnist
}


show.digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}