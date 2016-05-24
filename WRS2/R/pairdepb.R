pairdepb <- function(y, groups, blocks, tr = 0.2, nboot = 599){
#
  cols1 <- deparse(substitute(y))
  cols2 <- deparse(substitute(groups))
  cols3 <- deparse(substitute(blocks))
  dat <- data.frame(y, groups, blocks)
  colnames(dat) <- c(cols1, cols2, cols3)
  cl <- match.call()
  
  x <- reshape(dat, idvar = cols3, timevar = cols2, direction = "wide")[-1]  ## wide format
  grp <- c(1:length(x))
  
  alpha=.05
  grp=0
  if(is.data.frame(x)) x <- as.matrix(x)
  if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
  if(is.list(x)){
    if(sum(grp)==0)grp<-c(1:length(x))
    # put the data in an n by J matrix
    mat<-matrix(0,length(x[[1]]),length(grp))
    for (j in 1:length(grp))mat[,j]<-x[[grp[j]]]
  }
  if(is.matrix(x)){
    if(sum(grp)==0)grp<-c(1:ncol(x))
    mat<-x[,grp]
  }
  if(sum(is.na(mat)>=1))stop("Missing values are not allowed.")
  J<-ncol(mat)
  connum<-(J^2-J)/2
  bvec<-matrix(0,connum,nboot)
 # set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  #print("Taking bootstrap samples. Please wait.")
  data<-matrix(sample(nrow(mat),size=nrow(mat)*nboot,replace=TRUE),nrow=nboot)
  xcen<-matrix(0,nrow(mat),ncol(mat))
  for (j in 1:J)xcen[,j]<-mat[,j]-mean(mat[,j],tr) #Center data
  
  it<-0
  for (j in 1:J){
    for (k in 1:J){
      if(j<k){
        it<-it+1
        bvec[it,]<-apply(data,1,tsub,xcen[,j],xcen[,k],tr)
        # bvec is a connum by nboot matrix containing the bootstrap test statistics.
      }}}
  bvec<-abs(bvec)  #Doing two-sided confidence intervals
  icrit<-round((1-alpha)*nboot)
  critvec<-apply(bvec,2,max)
  critvec<-sort(critvec)
  crit<-critvec[icrit]
  psihat<-matrix(0,connum,5)
  dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper"))
  test<-matrix(NA,connum,4)
  dimnames(test)<-list(NULL,c("Group","Group","test","se"))
  it<-0
  for (j in 1:J){
    for (k in 1:J){
      if(j<k){
        it<-it+1
        estse<-yuend(mat[,j],mat[,k])$se
        dif<-mean(mat[,j],tr)-mean(mat[,k],tr)
        psihat[it,1]<-grp[j]
        psihat[it,2]<-grp[k]
        psihat[it,3]<-dif
        psihat[it,4]<-dif-crit*estse
        psihat[it,5]<-dif+crit*estse
        test[it,1]<-grp[j]
        test[it,2]<-grp[k]
        test[it,3]<-yuend(mat[,j],mat[,k])$test
        test[it,4]<-estse
      }}}
 
 fnames <- as.character(unique(groups))
 psihat1 <- cbind(psihat, test[,3], crit)
 
 result <- list(comp = psihat1, fnames = fnames, call = cl)
 class(result) <- "mcp1"
 result
}
