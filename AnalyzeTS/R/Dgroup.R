Dgroup <-
function(x,follow=NULL,r=2,answer=1,statistic="ALL"){

sosanh<-function(vector,x){
s=0
for(i in 1:length(vector)) if(vector[i]==x)s=s+1
if(s==0) kq<-"No" else kq<-"Yes"
kq
}

where<-function(vector,x){
id<-1:length(vector)
vt<-id[as.vector(vector)==x]
vt
}

is.wholenumber<-function(x,tol=.Machine$double.eps^0.5)
abs(x - round(x)) < tol

if(!is.wholenumber(r))stop("'r' must be a integer number!")
if(answer!=1 & answer!=2) stop("'answer' must be 1 or 2!")

if(is.null(follow)){
ans.statistic<-Descriptives(x,r=r,answer=answer,statistic=statistic)
}

if(!is.null(follow)){
if(!is.ts(x) & !is.vector(x)) stop("Sorry! We just statistic following group for 'time series' or 'vector' object.")

if(!is.list(follow)){
if(!is.factor(follow)) stop("'follow' must be a factor!")
if(length(follow)!=length(x)) stop("'x' and 'follow' are different about length!")
factor.fl<-"a"
for(i in 1:length(follow)) if(sosanh(factor.fl,as.vector(follow[i]))=="No") factor.fl<-cbind(factor.fl,as.vector(follow[i]))
factor.fl<-as.vector(factor.fl)[-1]
#statistic 1
temp<-c(1,2,3,4,5,2,3,6,7,2,5,6,7)
several<-Descriptives(temp,r=r,answer=1,statistic=statistic)
for(j in 1:length(factor.fl)){
locate<-where(follow,factor.fl[j])
several<-cbind(several,Descriptives(x[locate],r=r,answer=1,statistic=statistic))}
several<-several[,-1]
dimnames(several)[[2]]<-factor.fl
if(answer==2) several<-t(several)
ans.statistic<-several
#statistic 1 finish
}#follow is not a list


if(is.list(follow)){
if(length(follow)==1){
follow<-follow[[1]]
if(!is.factor(follow)) stop("'follow' must be a factor!")
if(length(follow)!=length(x)) stop("'x' and 'follow' are different about length!")
factor.fl<-"a"
for(i in 1:length(follow)) if(sosanh(factor.fl,as.vector(follow[i]))=="No") factor.fl<-cbind(factor.fl,as.vector(follow[i]))
factor.fl<-as.vector(factor.fl)[-1]
#statistic 1
temp<-c(1,2,3,4,5,2,3,6,7,2,5,6,7)
several<-Descriptives(temp,r=r,answer=1,statistic=statistic)
for(j in 1:length(factor.fl)){
locate<-where(follow,factor.fl[j])
several<-cbind(several,Descriptives(x[locate],r=r,answer=1,statistic=statistic))}
several<-several[,-1]
dimnames(several)[[2]]<-factor.fl
if(answer==2) several<-t(several)
ans.statistic<-several
#statistic 2 finish
}}#follow is a list 1 element

if(is.list(follow)){
if(length(follow)==2){
follow1<-follow[[1]]
follow2<-follow[[2]]
if(!is.factor(follow1)) stop("'follow[[1]]' must be a factor!")
if(!is.factor(follow2)) stop("'follow[[2]]' must be a factor!")
if(length(follow1)!=length(follow2)) stop("'follow[[1]]' and 'follow[[2]]' are different about length!")
if(length(follow1)!=length(x)) stop("'x' and 'follow[[1]]' are different about length!")
if(length(follow2)!=length(x)) stop("'x' and 'follow[[2]]' are different about length!")
dat<-data.frame(follow1,follow2,x)

levels1<-levels(follow1)

list.ans<-as.list(levels1)
names(list.ans)<-levels1

for(f1 in 1:length(levels1)){
locate1<-where(follow1,levels1[f1])
follow3<-follow2[locate1]
x3<-x[locate1]

factor.fl<-"a"
for(i in 1:length(follow3)) if(sosanh(factor.fl,as.vector(follow3[i]))=="No") factor.fl<-cbind(factor.fl,as.vector(follow3[i]))
factor.fl<-as.vector(factor.fl)[-1]

temp<-c(1,2,3,4,5,2,3,6,7,2,5,6,7)
several<-Descriptives(temp,r=r,answer=1,statistic=statistic)
for(j in 1:length(factor.fl)){
locate<-where(follow3,factor.fl[j])
several<-cbind(several,Descriptives(x3[locate],r=r,answer=1,statistic=statistic))}
several<-several[,-1]
dimnames(several)[[2]]<-factor.fl
list.ans[[f1]]<-several
}
if(answer==2){
for(i in 1:length(list.ans))
list.ans[[i]]<-t(list.ans[[i]])
}
ans.statistic<-list.ans
}}#follow is a list 2 element

if(is.list(follow))
if(length(follow)>2)stop("Sorry! This function just calculate statistic bigest 2 factor.")


}#finish follow != NULL
ans.statistic
}
