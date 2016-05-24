# stuart.maxwell.mh computes the marginal homogeneity test for
# a CxC matrix of assignments of objects to C categories or an
# nx2 or 2xn matrix of category scores for n data objects by two
# raters. The statistic is distributed as Chi-square with C-1
# degrees of freedom.

stuart.maxwell.mh<-function(x) {
 if(is.matrix(x) || is.data.frame(x)) {
  dimx<-dim(x)
  if(length(dimx) == 2) {
   if(dimx[1] != dimx[2]) {
    # if dimension lengths are unequal, assume it's a score matrix
    if(dimx[1] == 2) smx<-as.matrix(table(x[1,],x[2,]))
    # assume the matrix is nx2
    else smx<-as.matrix(table(x[,1],x[,2]))
   }
   else smx<-as.matrix(x)
   # get the marginals
   rowsums<-apply(smx,1,sum)
   colsums<-apply(smx,2,sum)
   equalsums<-rowsums == colsums
   if(any(equalsums)) {
    # dump any categories with perfect agreement
    smx<-smx[!equalsums,!equalsums]
    # bail out if too many categories have disappeared
    if(dim(smx)[1] < 2) stop("Too many equal marginals, cannot compute")
    # get new marginals
    rowsums<-apply(smx,1,sum)
    colsums<-apply(smx,2,sum)
   }
   # use K-1 marginals
   Kminus1<-length(rowsums)-1
   smd<-(rowsums-colsums)[1:Kminus1]
   smS<-matrix(0,nrow=Kminus1,ncol=Kminus1)
   for(i in 1:Kminus1) {
    for(j in 1:Kminus1) {
     if(i == j) smS[i,j]<-rowsums[i] + colsums[j] - 2 * smx[i,j]
     else smS[i,j]<--(smx[i,j]+smx[j,i])
    }
   }
   smstat<-t(smd)%*%solve(smS)%*%smd
   p=1-pchisq(smstat,Kminus1)
   stat.name<-paste("Chisq(",Kminus1,")",sep="")
   smmh<-structure(list(method="Stuart-Maxwell marginal homogeneity",
    subjects=sum(smx),raters=2,irr.name="Chisq",value=smstat,
    stat.name=stat.name,statistic=smstat,p.value=p),class="irrlist")
   class(smmh)<-"irrlist"
   return(smmh)
  }
  else cat("Dimension higher than 2, cannot compute\n")  
 }
 else {
  cat("Usage: stuart.maxwell.mh(x)\n")
  cat("\twhere x is an nx2 matrix or data frame of category scores for n objects\n")
  cat("\tor a CxC matrix or data frame of rater agreement on C categories\n")
 }
}

