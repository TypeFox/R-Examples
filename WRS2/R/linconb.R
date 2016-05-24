linconb<-function(x,con=0,tr=.2,alpha=.05,nboot=599,pr=TRUE,SEED=TRUE){
#
#   Compute a 1-alpha confidence interval for a set of d linear contrasts
#   involving trimmed means using the bootstrap-t bootstrap method.
#   Independent groups are assumed.
#
#   The data are assumed to be stored in x in list mode.  Thus,
#   x[[1]] contains the data for the first group, x[[2]] the data
#   for the second group, etc. Length(x)=the number of groups = J, say.
#
#   Missing values are automatically removed.
#
#   con is a J by d matrix containing the contrast coefficents of interest.
#   If unspecified, all pairwise comparisons are performed.
#   For example, con[,1]=c(1,1,-1,-1,0,0) and con[,2]=c(,1,-1,0,0,1,-1)
#   will test two contrasts: (1) the sum of the first two trimmed means is
#   equal to the sum of the second two, and (2) the difference between
#   the first two is equal to the difference between the trimmed means of
#   groups 5 and 6.
#
#   The default number of bootstrap samples is nboot=599
#
#   This function uses functions trimparts and trimpartt written for this
#   book.
#
#
#
#
if(is.data.frame(x))x=as.matrix(x)
#if(pr){
#print("Note: confidence intervals are adjusted to control FWE")
#print("But p-values are not adjusted to control FWE")
#}
con<-as.matrix(con)
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
J<-length(x)
for(j in 1:J){
xx<-x[[j]]
x[[j]]<-xx[!is.na(xx)] # Remove any missing values.
}
Jm<-J-1
d<-(J^2-J)/2
if(sum(con^2)==0){
con<-matrix(0,J,d)
id<-0
for (j in 1:Jm){
jp<-j+1
for (k in jp:J){
id<-id+1
con[j,id]<-1
con[k,id]<-0-1
}}}
if(nrow(con)!=length(x))stop("The number of groups does not match the number of contrast coefficients.")
bvec<-array(0,c(J,2,nboot))
if(SEED)set.seed(2) # set seed of random number generator so that
#             results can be duplicated.
#if(pr)print("Taking bootstrap samples. Please wait.")
nsam=matl(lapply(x,length))
for(j in 1:J){
paste("Working on group ",j)
xcen<-x[[j]]-mean(x[[j]],tr)
data<-matrix(sample(xcen,size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
bvec[j,,]<-apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
#                     contains the bootstrap trimmed means, the second row
#                     contains the bootstrap squared standard errors.
}
m1<-bvec[,1,]  # J by nboot matrix containing the bootstrap trimmed means
m2<-bvec[,2,]  # J by nboot matrix containing the bootstrap sq. se.
boot<-matrix(0,ncol(con),nboot)
for (d in 1:ncol(con)){
top<-apply(m1,2,trimpartt,con[,d])
#            A vector of length nboot containing psi hat values
consq<-con[,d]^2
bot<-apply(m2,2,trimpartt,consq)
boot[d,]<-abs(top)/sqrt(bot)
}
testb<-apply(boot,2,max)
ic<-floor((1-alpha)*nboot)
testb<-sort(testb)
psihat<-matrix(0,ncol(con),4)
test<-matrix(0,ncol(con),4)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper"))
dimnames(test)<-list(NULL,c("con.num","test","se","p.value"))
for (d in 1:ncol(con)){
test[d,1]<-d
psihat[d,1]<-d
testit<-lincon1(x,con[,d],tr,pr=FALSE)
test[d,2]<-testit$test[1,2]
pval<-mean((abs(testit$test[1,2])<boot[d,]))
test[d,4]<-pval
psihat[d,3]<-testit$psihat[1,2]-testb[ic]*testit$test[1,4]
psihat[d,4]<-testit$psihat[1,2]+testb[ic]*testit$test[1,4]
psihat[d,2]<-testit$psihat[1,2]
test[d,3]<-testit$test[1,4]
}
list(n=nsam,psihat=psihat,test=test,crit=testb[ic],con=con)
}
