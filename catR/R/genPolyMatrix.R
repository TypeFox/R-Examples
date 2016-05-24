genPolyMatrix<-function(items=100, nrCat=3, model="GRM",seed=1,same.nrCat=FALSE,cbControl=NULL){
if (sum(model==c("GRM","MGRM","PCM","GPCM","RSM","NRM"))==0) stop("invalid 'model' name'",call.=FALSE)
set.seed(seed)
if (same.nrCat | model=="MGRM" | model=="RSM") gj<-rep(nrCat-1,items)
else {
gj<-rpois(items,nrCat-1)
gj[gj>nrCat-1]<-nrCat-1
gj[gj<1]<-1
}
if (model=="GRM" | model=="GPCM"){
res<-matrix(NA,items,(max(gj)+1))
alphaj<-rlnorm(items,0,0.1225)
for (i in 1:items){
pars<-rnorm(gj[i])
if (model=="GRM") pars<-sort(pars)
res[i,1:(length(pars)+1)]<-c(alphaj[i],pars)
}
name<-"alphaj"
for (i in 1:max(gj)) {
if (model=="GRM") name<-c(name,paste("betaj",i,sep=""))
else name<-c(name,paste("deltaj",i,sep=""))
}
colnames(res)<-name
}

if (model=="MGRM"){
res<-matrix(NA,items,(max(gj)+2))
alphaj<-rlnorm(items,0,0.1225)
bj<-rnorm(items)
pars<-sort(rnorm(max(gj)),decreasing=TRUE)
for (i in 1:nrow(res)) res[i,]<-c(alphaj[i],bj[i],pars)
name<-c("alphaj","bj")
for (i in 1:max(gj)) name<-c(name,paste("c",i,sep=""))
colnames(res)<-name
}

if (model=="PCM"){
res<-matrix(NA,items,max(gj))
for (i in 1:items){
pars<-rnorm(gj[i])
res[i,1:length(pars)]<-pars
}
name<-NULL
for (i in 1:max(gj)) name<-c(name,paste("deltaj",i,sep=""))
colnames(res)<-name
}

if (model=="RSM"){
res<-matrix(NA,items,(max(gj)+1))
lambdaj<-rnorm(items)
deltaj<-rnorm(max(gj))
for (i in 1:items) res[i,]<-c(lambdaj[i],deltaj)
name<-c("lambdaj")
for (i in 1:max(gj)) name<-c(name,paste("delta",i,sep=""))
colnames(res)<-name
}

if (model=="NRM"){
res<-matrix(NA,items,2*max(gj))
for (i in 1:items){
alphaj<-rlnorm(gj[i],0,0.1225)
cj<-rnorm(gj[i])
v<-NULL
for (t in 1:gj[i]) v<-c(v,alphaj[t],cj[t])
res[i,1:(2*gj[i])]<-v
}
name<-NULL
for (i in 1:max(gj)) name<-c(name,paste("alpha",i,sep=""),paste("c",i,sep=""))
colnames(res)<-name
}
if (!is.null(cbControl)){
pr<-rmultinom(items,1,cbControl$props)
Group<-NULL
for (i in 1:items) Group[i]<-cbControl$names[pr[,i]==1]
res<-data.frame(round(res,3),Group)
#res<-res[order(res[,ncol(res)]),]
}
else res<-data.frame(round(res,3))
return(res)
}


