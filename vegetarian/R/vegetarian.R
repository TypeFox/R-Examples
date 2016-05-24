####Vegetarian Library
#______________________________________________________
normalize.rows<-function(abundance.matrix){
#Makes it so each row sums to one
if (is.data.frame(abundance.matrix)) {abundance.matrix <- as.matrix(abundance.matrix)}


if(is.matrix(abundance.matrix)){
norm.abundance<-matrix(NA,nrow=dim(abundance.matrix)[1],ncol=dim(abundance.matrix)[2])
for(i in 1:dim(abundance.matrix)[1]){
sum.row<-sum(abundance.matrix[i,])
if(sum.row!=0){norm.abundance[i,]<-abundance.matrix[i,]/sum.row}
if(sum.row==0){norm.abundance[i,]<-abundance.matrix[i,]}
}
}

if(is.vector(abundance.matrix)){
sum.row<-sum(abundance.matrix)
if(sum.row!=0){norm.abundance<-abundance.matrix/sum.row}
if(sum.row==0){norm.abundance<-abundance.matrix}
}
norm.abundance
}

 
#______________________________________________________

p.q.sum<-function(p,q=1){
temp.sum<-0
for(i in 1:length(p)){
if(p[i]>0){
if(q!=1){temp.sum<-temp.sum+p[i]^q}
if(q==1){temp.sum<-temp.sum+p[i]*log(p[i])}
#if q=1, then we are doing Shannon type diversities, e.g. Jost 2008 equation 17a
}
}
temp.sum
}

 
#______________________________________________________
bootstrap<-function(abundances,s.sizes=NULL,num.iter=100,func,func.arg="blank",sim.pop=FALSE,sim.par=FALSE){
#func is the function (e.g. d)
#arg is a list of the arguments the function needs (e.g. list(q=1,wts=FALSE) )
require(stats)
num.sites<-dim(abundances)[1]
num.spec<-dim(abundances)[2]
if(is.null(s.sizes)){s.sizes<-apply(abundances,1,sum)} 
#if s.sizes is not specified, then assumes scores are counts, and sum or rows is sample sizes
boots<-array(NA,dim=c(num.sites,num.spec,num.iter))
simulated.param<-vector(length=num.iter)


for(j in 1:num.sites){
temp.boot<-rmultinom(n=num.iter, size=s.sizes[j],prob=abundances[j,])
boots[j,,]<-temp.boot
}

for(i in 1:num.iter){
if(func.arg[1]=="blank"){simulated.param[i]<-do.call(what=func,arg=list(boots[,,i]))}
if(func.arg[1]!="blank"){simulated.param[i]<-do.call(what=func,arg=c(abundances=list(boots[,,i]),func.arg))}
}


variance<-( sum( (simulated.param-mean(simulated.param) )^2) )/length(simulated.param)
stderr<-sqrt(variance)

if(sim.pop==FALSE & sim.par==FALSE){
OUT<-list(StdErr=stderr)
}
if(sim.pop==FALSE & sim.par==TRUE){
OUT<-list(Simulated.Parameters=simulated.param,StdErr=stderr)
}

if(sim.pop==TRUE & sim.par==FALSE){
OUT<-list(Simulated.Populations=boots,StdErr=stderr)
}

if(sim.pop==TRUE & sim.par==TRUE){
OUT<-list(Simulated.Populations=boots,Simulated.Parameters=simulated.param,StdErr=stderr)
}

OUT
}

#____________________________________________________________
d<-function(abundances,lev="alpha",wts=FALSE,q=1,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
#boot.arg is a list of arguments to feed bootstrap (e.g. boot.arg=list(s.sizes=c(10,12,10,20,6),num.iter=50) )
if(is.data.frame(abundances)){abundances<-as.matrix(abundances)}

if(wts[1]==FALSE){
if(is.matrix(abundances)){wts<-rep(1,dim(abundances)[1])}
if(is.vector(abundances)){wts<-1}
}

if(lev=="alpha"){
if(q=="shannon"){q<-1}
p<-normalize.rows(abundances)
w<-normalize.rows(wts)

if(q!=1){
#Jost 2008 equation 11a
numerator<-0
if(is.matrix(abundances)){
for(j in 1:dim(p)[1]){
numerator<-numerator+(w[j]^q)*(p.q.sum(p[j,],q=q))
}
}

if(is.vector(abundances)){
numerator<-p.q.sum(p,q=q)
}

D.VALUE<- ( ( numerator/(sum(w^q)) )^(1/(1-q)) )
}

if(q==1){
#Shannon Case; Jost 2008 equation 11b

if(is.matrix(abundances)){
numerator<-0
for(j in 1:dim(p)[1]){
numerator<-numerator+((-1)*w[j])*(p.q.sum(p[j,],q=q))
#wts (w[j]) should sum to one, so this essentially is a mean alpha over all rows
}
}

if(is.vector(abundances)){
numerator<-(-1)*p.q.sum(p,q=q)
}
D.VALUE<-exp(numerator)
}
}

if(lev=="beta"){D.VALUE<-d(abundances,lev="gamma",wts=wts,q=q)/d(abundances,lev="alpha",wts=wts,q=q)}

if(lev=="gamma"){
#gamma diveristy is just alpha metric applied to pooled sites
if(wts[1]==FALSE){
if(is.matrix(abundances)){wts<-rep(1,dim(abundances)[1])}
if(is.vector(abundances)){wts<-1}
}
w<-normalize.rows(wts)

D.VALUE<-d(apply(normalize.rows(abundances)*w,2,sum),lev="alpha",q=q)
}



if(boot==TRUE){
func.arg<-list(wts=wts,q=q,lev=lev)
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(abundances),func=d,func.arg=list(func.arg),boot.arg))
output<-list(D.Value=D.VALUE,StdErr=boot.out$StdErr)
}

if(boot==FALSE){output<-D.VALUE}

output

}

#__________________________________________________________________________
H<-function(abundances,lev="alpha",wts=FALSE,q=1,HCDT=FALSE,gini=FALSE,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
#based on Jost 2008 table 1
if(is.data.frame(abundances)){abundances<-as.matrix(abundances)}


#richness
if(q==0){(H.VALUE<-d(abundances,lev=lev,wts=wts,q=q))}

#Shannon entropy
if(q==1|q=="shannon"){H.VALUE<-log(d(abundances,lev=lev,wts=wts,q=q))}

#Simpson concentration
if(q==2){H.VALUE<-(1/(d(abundances,lev=lev,wts=wts,q=q)))
if(gini==TRUE){H.VALUE<-(1-H.VALUE)}
}


if(q!=0&q!=1&q!=2&HCDT==FALSE){

#Renyi entropy
if(HCDT==FALSE){H.VALUE<-log(d(abundances,lev=lev,wts=wts,q=q))}

#HCDT entropy
if(HCDT==TRUE){
Da<-d(abundances,lev=lev,wts=wts,q=q)
H.VALUE<-( (Da^(1-q) - 1)/(1-q) )
}
}

if(boot==TRUE){
func.arg<-list(wts=wts,q=q,lev=lev)
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(abundances),func=H,func.arg=list(func.arg),boot.arg))
output<-list(H.Value=H.VALUE,StdErr=boot.out$StdErr)
}

if(boot==FALSE){output<-H.VALUE}

output

}

#________________________________________________________________________

 
M.homog<-function(abundances,abundances2=NULL,q=1,std=FALSE,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
#can be given as two vectors, or a single matrix with two rows
if(is.data.frame(abundances)){abundances<-as.matrix(abundances)}
if(is.data.frame(abundances2)){abundances2<-as.vector(abundances2)}

if(!is.null(abundances2)){
abundances<-rbind(abundances,abundances2)
}
N<-dim(abundances)[1]

#Jost 2008 equation 19
M<-1/d(abundances,lev="beta",q=q)

#Jost 2008 equation 20
if(std==TRUE){M<-(M-(1/N))/(1-(1/N))}

if(boot==TRUE){
func.arg<-list(q=q)
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(abundances),func=M.homog,func.arg=list(func.arg),boot.arg))
output<-list(M=M,StdErr=boot.out$StdErr)
}
if(boot==FALSE){output<-M}


output
}
#______________________________________________________

Rel.homog<-function(abundances,abundances2=NULL,wts=FALSE,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
if(!is.null(abundances2)){
abundances<-rbind(abundances,abundances2)
}

if(wts[1]==FALSE){wts<-rep(1,dim(abundances)[1])}

w<-normalize.rows(wts)

#Jost 2008 equation 21
temp.sum<-0
j<-1

for(j in 1:dim(abundances)[1]){
temp.sum<-(temp.sum+w[j]*log(w[j]))
}
Dw<-exp(-temp.sum)
Db<-d(abundances,lev="beta",q=1)
HOMOG<-((1/Db) - (1/Dw))/(1 - (1/Dw))

if(boot==TRUE){
func.arg<-list(wts=wts)
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(abundances),func=Rel.homog,func.arg=list(func.arg),boot.arg))
output<-list(REL.HOMOG=HOMOG,StdErr=boot.out$StdErr)
}

if(boot==FALSE){output<-HOMOG}
output
}

#______________________________________________________
similarity<-function(abundances,abundances2=NULL,q=1,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
if(!is.null(abundances2)){
abundances<-rbind(abundances,abundances2)
}

N<-dim(abundances)[1]
if(q=="shannon"){q<-1}

#Jost 2008 equation 24
#If q==1, this is only valid for 2 equally weighted communities
if(q==1){S<-(log(N)-H(abundances,lev="beta",wts=FALSE,q=1))/log(N)}

#Jost 2008 equation 23
if(q!=1){
N<-dim(abundances)[1]
Db<-d(abundances,lev="beta",q=q)
S<-(((1/Db)^(q-1) - (1/N)^(q-1))/(1 - (1/N)^(q-1)))
}

if(boot==TRUE){
func.arg<-list(q=q)
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(abundances),func=similarity,func.arg=list(func.arg),boot.arg))
output<-list(Simlarity=S,StdErr=boot.out$StdErr)
}

if(boot==FALSE){output<-S}

output
}

#______________________________________________________


turnover<-function(abundances,abundances2=NULL,q=1,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
if(is.vector(abundances)&!is.null(abundances2)){
abundances<-rbind(abundances,abundances2)
}

Db<-d(abundances,lev="beta",q=q)
N<-dim(abundances)[1]

#Jost 2008 equation 25
Tv<-(Db-1)/(N-1)

if(boot==TRUE){
func.arg<-list(q=q)
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(abundances),func=turnover,func.arg=list(func.arg),boot.arg))
output<-list(Turnover=Tv,StdErr=boot.out$StdErr)
}

if(boot==FALSE){output<-Tv}

output
}

#______________________________________________________

sim.table<-function(abundances,q=1,labels=FALSE,half=TRUE,diag=TRUE,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
num.sites<-dim(abundances)[1]
ifelse(diag==TRUE,k<-0,k<-1)
func.arg<-list(q=q)

if(length(labels)>1){
#assumes site names handed as vector, and no labels in main matrix
names<-list(labels,labels)
scores<-abundances
}

if(length(labels)==1){

if(labels==TRUE){
#assumes site names are the first column of the data matrix
names<-list(abundances[,1],abundances[,1])
scores<-as.matrix(abundances[,-1])
}

if(labels==FALSE){
#assumes no site names given anywhere
names<-list(seq(1:num.sites),seq(1:num.sites))
scores<-as.matrix(abundances)
}
}

sim.matrix<-matrix(NA,nrow=num.sites,ncol=num.sites,dimnames=names)

if(half==TRUE){
for(i in (1+k):num.sites){
for(j in 1:(i-k)){

if(boot==FALSE){sim.matrix[i,j]<-similarity(scores[i,],scores[j,],q=q)}
if(boot==TRUE){
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(rbind(scores[i,],scores[j,])),func=similarity,func.arg=list(func.arg),boot.arg))
sim.matrix[i,j]<-boot.out$StdErr
}

}
}
}

if(half==FALSE){
for(i in 1:num.sites){
for(j in 1:num.sites){

if(boot==FALSE){sim.matrix[i,j]<-similarity(scores[i,],scores[j,],q=q)}
if(boot==TRUE){
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(rbind(scores[i,],scores[j,])),func=similarity,func.arg=list(func.arg),boot.arg))
sim.matrix[i,j]<-boot.out$StdErr
}

}
}
}
sim.matrix
}

#______________________________________________________


sim.groups<-function(abundances1,abundances2,q=1,labels=FALSE,boot=FALSE,boot.arg=list(s.sizes=NULL,num.iter=100)){
func.arg<-list(q=q)

num.sites1<-dim(abundances1)[1]

if(labels==TRUE){
abundance.1<-as.matrix(abundances1[,-1])
abundance.2<-as.matrix(abundances2[,-1])
scores<-rbind(abundance.1,abundance.2)
}

if(labels==FALSE){
abundance.1<-as.matrix(abundances1)
abundance.2<-as.matrix(abundances2)
scores<-rbind(abundance.1,abundance.2)
}


num.sites.tot<-dim(scores)[1]
within.1<-NA
within.2<-NA
btwn<-NA


for(i in 2:num.sites.tot){
for(j in 1:(i-1)){

if(boot==FALSE){temp.sim<-similarity(scores[i,],scores[j,],q=q)}
if(boot==TRUE){
boot.out<-do.call(what=bootstrap,arg=c(abundances=list(rbind(scores[i,],scores[j,])),func=similarity,func.arg=list(func.arg),boot.arg))
temp.sim<-boot.out$StdErr
}

if(i<=num.sites1){within.1<-c(within.1,temp.sim)}

if(i>num.sites1){
if(j>num.sites1){within.2<-c(within.2,temp.sim)}
if(j<=num.sites1){btwn<-c(btwn,temp.sim)}
}
rm(temp.sim)
}
}
within.1<-within.1[-1]
within.2<-within.2[-1]
btwn<-btwn[-1]

list(within.1=within.1,within.2=within.2,between=btwn)
} 

