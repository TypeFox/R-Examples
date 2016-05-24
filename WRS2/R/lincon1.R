lincon1<-function(x,con=0,tr=.2,alpha=.05,pr=TRUE,crit=NA,SEED=TRUE,KB=FALSE){
#
#  A heteroscedastic test of d linear contrasts using trimmed means.
#
#  The data are assumed to be stored in $x$ in list mode.
#  Length(x) is assumed to correspond to the total number of groups.
#  It is assumed all groups are independent.
#
#  con is a J by d matrix containing the contrast coefficients that are used.
#  If con is not specified, all pairwise comparisons are made.
#
#  Missing values are automatically removed.
#
#  To apply the Kaiser-Bowden method, use the function kbcon
#
if(tr==.5)stop("Use the R function medpb to compare medians")
if(is.data.frame(x))x=as.matrix(x)
if(KB)stop("Use the function kbcon")
flag<-T
if(alpha!= .05 && alpha!=.01)flag<-F
if(is.matrix(x))x<-listm(x)
if(!is.list(x))stop("Data must be stored in a matrix or in list mode.")
con<-as.matrix(con)
J<-length(x)
sam=NA
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
sam[j]=length(x[[j]])
h[j]<-length(x[[j]])-2*floor(tr*length(x[[j]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-((length(x[[j]])-1)*winvar(x[[j]],tr))/(h[j]*(h[j]-1))
xbar[j]<-mean(x[[j]],tr)
}
if(sum(con^2)==0){
CC<-(J^2-J)/2
#if(CC>28)print("For faster execution time but less power, use kbcon")
psihat<-matrix(0,CC,6)
dimnames(psihat)<-list(NULL,c("Group","Group","psihat","ci.lower","ci.upper",
"p.value"))
test<-matrix(NA,CC,6)
dimnames(test)<-list(NULL,c("Group","Group","test","crit","se","df"))
jcom<-0
for (j in 1:J){
for (k in 1:J){
if (j < k){
jcom<-jcom+1
test[jcom,3]<-abs(xbar[j]-xbar[k])/sqrt(w[j]+w[k])
sejk<-sqrt(w[j]+w[k])
test[jcom,5]<-sejk
psihat[jcom,1]<-j
psihat[jcom,2]<-k
test[jcom,1]<-j
test[jcom,2]<-k
psihat[jcom,3]<-(xbar[j]-xbar[k])
df<-(w[j]+w[k])^2/(w[j]^2/(h[j]-1)+w[k]^2/(h[k]-1))
test[jcom,6]<-df
psihat[jcom,6]<-2*(1-pt(test[jcom,3],df))
if(!KB){
if(CC>28)flag=F
if(flag){
if(alpha==.05)crit<-smmcrit(df,CC)
if(alpha==.01)crit<-smmcrit01(df,CC)
}
if(!flag || CC>28)crit<-smmvalv2(dfvec=rep(df,CC),alpha=alpha,SEED=SEED)
}
if(KB)crit<-sqrt((J-1)*(1+(J-2)/df)*qf(1-alpha,J-1,df))
test[jcom,4]<-crit
psihat[jcom,4]<-(xbar[j]-xbar[k])-crit*sejk
psihat[jcom,5]<-(xbar[j]-xbar[k])+crit*sejk
}}}}
if(sum(con^2)>0){
if(nrow(con)!=length(x)){
stop("The number of groups does not match the number of contrast coefficients.")
}
psihat<-matrix(0,ncol(con),5)
dimnames(psihat)<-list(NULL,c("con.num","psihat","ci.lower","ci.upper",
"p.value"))
test<-matrix(0,ncol(con),5)
dimnames(test)<-list(NULL,c("con.num","test","crit","se","df"))
df<-0
for (d in 1:ncol(con)){
psihat[d,1]<-d
psihat[d,2]<-sum(con[,d]*xbar)
sejk<-sqrt(sum(con[,d]^2*w))
test[d,1]<-d
test[d,2]<-sum(con[,d]*xbar)/sejk
df<-(sum(con[,d]^2*w))^2/sum(con[,d]^4*w^2/(h-1))
if(flag){
if(alpha==.05)crit<-smmcrit(df,ncol(con))
if(alpha==.01)crit<-smmcrit01(df,ncol(con))
}
if(!flag)crit<-smmvalv2(dfvec=rep(df,ncol(con)),alpha=alpha,SEED=SEED)
test[d,3]<-crit
test[d,4]<-sejk
test[d,5]<-df
psihat[d,3]<-psihat[d,2]-crit*sejk
psihat[d,4]<-psihat[d,2]+crit*sejk
psihat[d,5]<-2*(1-pt(abs(test[d,2]),df))
}
}

list(n=sam,test=test,psihat=psihat)
}
