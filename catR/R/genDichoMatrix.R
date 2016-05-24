genDichoMatrix<-function(items=100, cbControl=NULL, model="4PL", aPrior=c("norm",1,0.2),bPrior=c("norm",0,1),cPrior=c("unif",0,0.25),dPrior=c("unif",0.75,1), seed=1){
testB<-switch(bPrior[1],norm=1,unif=2)
if (is.null(testB)) stop("Prior distribution for item difficulties must be either 'norm' or 'unif'",call.=FALSE)
testA<-switch(aPrior[1],norm=1,lnorm=2,unif=3)
if (is.null(testA)) stop("Prior distribution for item discriminations must be either 'norm', 'lnorm' or 'unif'",call.=FALSE)
testC<-switch(cPrior[1],beta=1,unif=2)
if (is.null(testC)) stop("Prior distribution for item lower asymptotes must be either 'beta' or 'unif'",call.=FALSE)
testD<-switch(dPrior[1],beta=1,unif=2)
if (is.null(testD)) stop("Prior distribution for item upper asymptotes must be either 'beta' or 'unif'",call.=FALSE)
set.seed(seed)
b<-switch(bPrior[1],norm=rnorm(items,as.numeric(bPrior[2]),as.numeric(bPrior[3])),unif=runif(items,as.numeric(bPrior[2]),as.numeric(bPrior[3])))
if (model!="1PL") a<-switch(aPrior[1],norm=rnorm(items,as.numeric(aPrior[2]),as.numeric(aPrior[3])),lnorm=rlnorm(items,as.numeric(aPrior[2]),as.numeric(aPrior[3])),unif=runif(items,as.numeric(aPrior[2]),as.numeric(aPrior[3])))
else a<-rep(1,items) 
if (model=="3PL" | model=="4PL") c<-switch(cPrior[1],beta=rbeta(items,as.numeric(cPrior[2]),as.numeric(cPrior[3])),unif=runif(items,as.numeric(cPrior[2]),as.numeric(cPrior[3])))
else c<-rep(0,items)
if (model=="4PL") d<-switch(dPrior[1],beta=rbeta(items,as.numeric(dPrior[2]),as.numeric(dPrior[3])),unif=runif(items,as.numeric(dPrior[2]),as.numeric(dPrior[3])))
else d<-rep(1,items)
if (!is.null(cbControl)){
props<-cbControl$props
names<-cbControl$names
props<-props/sum(props)
nr<-round(items*props)
f<-function(items,nr){
res<-matrix(0,length(nr),items)
for (i in 1:items){
res[,i]<-rmultinom(1,1,nr)
nr<-nr-res[,i]
}
return(res)}
pr<-f(items,nr)
g<-function(x) (1:length(x))[x==1]
it<-apply(pr,2,g)
groups<-NULL
for (i in 1:items) groups[i]<-names[it[i]]
groups<-as.factor(groups)
RES<-data.frame(a,b,c,d,Group=groups)
}
else RES<-data.frame(a,b,c,d)
return(RES)}








