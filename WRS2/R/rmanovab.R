rmanovab <- function(y, groups, blocks, tr = 0.2, nboot = 599){
  
  cols1 <- deparse(substitute(y))
  cols2 <- deparse(substitute(groups))
  cols3 <- deparse(substitute(blocks))
  dat <- data.frame(y, groups, blocks)
  colnames(dat) <- c(cols1, cols2, cols3)
  cl <- match.call()
  
  x <- reshape(dat, idvar = cols3, timevar = cols2, direction = "wide")[-1]  ## wide format
  
  alpha <- 0.05
  grp <- 0
  if(is.data.frame(x)) x=as.matrix(x)
  
  
  if(is.matrix(x)){
    if(sum(grp)==0)grp<-c(1:ncol(x))
    mat<-x[,grp]
  }
  mat=elimna(mat)
  J<-ncol(mat)
  connum<-(J^2-J)/2
  bvec<-matrix(0,connum,nboot)
 # if (SEED) set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  data<-matrix(sample(nrow(mat),size=nrow(mat)*nboot,replace=TRUE),nrow=nboot)
  xcen<-matrix(0,nrow(mat),ncol(mat))
  for (j in 1:J)xcen[,j]<-mat[,j]-mean(mat[,j],tr) #Center data
  bvec<-apply(data,1,tsubrmanovab,xcen,tr)
  # bvec is vector of nboot  bootstrap test statistics.
  icrit<-round((1-alpha)*nboot)
  bvec<-sort(bvec)
  crit<-bvec[icrit]
  test<-rmanovatemp(mat,tr,grp)$test
  result <- list(test = test, crit = crit, call = cl)
  class(result) <- c("rmanovab")
  result
}
