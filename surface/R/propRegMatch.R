propRegMatch <-
function(fit1,fit2,internal=FALSE){

if(class(fit1)=="list")fit1<-fit1[[1]]
if(class(fit2)=="list")fit2<-fit2[[1]]

if(internal){ 
	fit1<-as(fit1,"data.frame")[,5,drop=F];fit2<-as(fit2,"data.frame")[,5,drop=F]
	n<-dim(fit1)[1]
}else{
	if(class(fit1)%in%c("character","factor"))fit1<-data.frame(regs=fit1)
	if(class(fit1)=="hansentree"){
		fit1<-as(fit1,"data.frame")
		tips<-((dim(fit1)[1]+1)/2):dim(fit1)[1];n<-length(tips)
		fit1<-data.frame(regs=fit1[tips,5],row.names=fit1$labels[tips])
		}
	if(class(fit2)%in%c("character","factor"))fit2<-data.frame(regs=fit2)
	if(class(fit2)=="hansentree"){
		fit2<-as(fit2,"data.frame")
		tips<-((dim(fit2)[1]+1)/2):dim(fit2)[1];n<-length(tips)
		fit2<-data.frame(regs=fit2[tips,5],row.names=fit2$labels[tips])
		}

	if(any(rownames(fit1)%in%rownames(fit2)==F)|any(rownames(fit2)%in%rownames(fit1)==F))stop("`fit1` and `fit2` must have identical names/labels")
	fit2<-fit2[rownames(fit1),,drop=F]	
}

mx1<-mx2<-matrix(NA,n,n,dimnames=list(rownames(fit1),rownames(fit1)))
for(j in 1:n){
	for(k in 1:n){
		if(j>k){
			mx1[j,k]<-as.numeric(fit1[j,1]==fit1[k,1])
			mx2[j,k]<-as.numeric(fit2[j,1]==fit2[k,1])
		}}}

pmatch<-sum(mx1==mx2,na.rm=T)/(n*(n-1)/2)
pmatch
}
