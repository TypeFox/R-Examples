CheckEstimateTFDRFormat <- function(dat,type){
  if(type!="multiclass" && type!="survival" && type!="regression" && type!="two class") stop("type must be one of 'multiclass', 'regression', 'two class', or 'survival'")
  if(is.null(dat$x)||!is.matrix(dat$x)) stop("dat$x must be a matrix....")
  p <- nrow(dat$x)
  n <- ncol(dat$x)
  if(length(dat$y)!=n) stop("dat$y must be a vector of length equal to the # of columns of dat$x....")
  if((sum(is.na(dat$x))+sum(is.na(dat$y)))>0) stop("Error... one of your inputs contains NA's....")
  if(type=="survival") if(length(dat$censoring.status)!=n || sum(dat$censoring.status!=1 & dat$censoring.status!=0)>0) stop("dat$censoring.status must be a vector of 1's and 0's of length equal to the # of columns of dat$x....")
  if(type=="survival") if(sum(dat$y<0)>0) stop("Survival times must be positive....")
  if(type=="two class") if(sum(dat$y!=1 & dat$y!=2)>0) stop("For two-class outcome, dat$y must take values 1 or 2....")
   if(type=="multiclass") if(sum(sort(unique(dat$y))!=(1:length(unique(dat$y))))>0 || length(unique(dat$y))<2){
     stop("For multiclass outcome, dat$y must take on values 1,2,....")
   }
} 
