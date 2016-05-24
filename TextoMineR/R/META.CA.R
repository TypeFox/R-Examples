META.CA <-
function(x,naxes=5,axe.x=1,axe.y=2,lmd=3,lmk=3, main="Metakeys & Metadocs", graph = TRUE){

Ccontr<-x$col$contrib
Ccoor<-x$col$coord
Rcontr<-x$row$contrib
Rcoor<-x$row$coord
eigen<-x$eig

if(naxes>=min(nrow(Rcontr),nrow(Ccontr)))
naxes=min((nrow(Rcontr)-1),(nrow(Ccontr)-1))

# Computing the metakeys: words with contribuions over lmk*average_contribution of words  
Metakeys<-vector(mode="list",length=naxes)
for (i in 1:naxes){
	Metakeys[[i]]=vector(mode="list",length=2)}
for(i in 1:naxes){
	Metakeys[[i]][[1]]<-sort(Ccontr[which(Ccontr[,i]>lmk*mean(Ccontr[,i])&Ccoor[,i]>0),i],decreasing=TRUE)
	if (length(Metakeys[[i]][[1]])==1) names(Metakeys[[i]][[1]])<-rownames(Ccontr)[which(Ccontr[,i]%in%sort(Ccontr[which(Ccontr[,i]>lmk*mean(Ccontr[,i])&Ccoor[,i]>0),i],decreasing=TRUE))]
	Metakeys[[i]][[2]]<-sort(Ccontr[which(Ccontr[,i]>lmk*mean(Ccontr[,i])&Ccoor[,i]<0),i],decreasing=TRUE)
	if (length(Metakeys[[i]][[2]])==1) names(Metakeys[[i]][[2]])<-rownames(Ccontr)[which(Ccontr[,i]%in%sort(Ccontr[which(Ccontr[,i]>lmk*mean(Ccontr[,i])&Ccoor[,i]<0),i],decreasing=TRUE))]
}

# Computing the metadocs : documents/answers with contribuions over lmd*average_contribution of documents 
Metadocs<-vector(mode="list",length=naxes)
for (i in 1:naxes){
	Metadocs[[i]]=vector(mode="list",length=2)}
for(i in 1:naxes){
	Metadocs[[i]][[1]]<-sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]>0),i],decreasing=TRUE)
	if (length(Metadocs[[i]][[1]])==1) names(Metadocs[[i]][[1]])<-rownames(Rcontr)[which(Rcontr[,i]%in%sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]>0),i],decreasing=TRUE))]
	Metadocs[[i]][[2]]<-sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]<0),i],decreasing=TRUE)
	if (length(Metadocs[[i]][[2]])==1) names(Metadocs[[i]][[2]])<-rownames(Rcontr)[which(Rcontr[,i]%in%sort(Rcontr[which(Rcontr[,i]>lmd*mean(Rcontr[,i])&Rcoor[,i]<0),i],decreasing=TRUE))]
}
# metakeys and metadocs
	 Metakeys.MetaDocs<-vector(mode="list",naxes)
       names( Metakeys.MetaDocs)<-paste("DIM",1:naxes,sep="")
	for (i in 1:naxes){
		X <-list(names(Metakeys[[i]][[1]]),names(Metadocs[[i]][[1]]),
                      names(Metakeys[[i]][[2]]),names(Metadocs[[i]][[2]]))
            names(X)<-c("Metakeys+","Metadocs+","Metakeys-","Metadocs-")
           Metakeys.MetaDocs[[i]]<-X 
          }

 #Dimension of the words
    mkeys<-numeric()
	for (i in 1:naxes){
		for (j in 1:2){
			mkeys<-c(mkeys,names(Metakeys[[i]][[j]]))
		}
	}
	mkeys<-as.factor(mkeys)
	levels(mkeys)
	matmkeys<-t(rbind(summary(mkeys,maxsum=10000),rep(naxes,length(levels(mkeys))),round(summary(mkeys,maxsum=10000)*100/naxes,2)))
	colnames(matmkeys)<-c("Dim","Total.Dim","%Dim")
	DimensionWord<-matmkeys[order(matmkeys[,1],decreasing = T),]
if(graph){
# Graph of metakeys and metadocs on the "axes" of the CA
axes<-c(axe.x,axe.y)

# Identify metakeys of the dimensions "axes"
ClusPal1<-which(Ccontr[,axes[1]]%in%Metakeys[[axes[1]]][[1]])
ClusPal2<-which(Ccontr[,axes[1]]%in%Metakeys[[axes[1]]][[2]])
ClusPal3<-which(Ccontr[,axes[2]]%in%Metakeys[[axes[2]]][[1]])
ClusPal4<-which(Ccontr[,axes[2]]%in%Metakeys[[axes[2]]][[2]])
ClusPal<-c(ClusPal1,ClusPal2,ClusPal3,ClusPal4)
ClusPal<-ClusPal[!duplicated(ClusPal)]
	
# Identify metadocs of the dimensions "axes"
ClusDoc1<-which(Rcontr[,axes[1]]%in%Metadocs[[axes[1]]][[1]])
ClusDoc2<-which(Rcontr[,axes[1]]%in%Metadocs[[axes[1]]][[2]])
ClusDoc3<-which(Rcontr[,axes[2]]%in%Metadocs[[axes[2]]][[1]])
ClusDoc4<-which(Rcontr[,axes[2]]%in%Metadocs[[axes[2]]][[2]])
ClusDoc<-c(ClusDoc1,ClusDoc2,ClusDoc3,ClusDoc4)
ClusDoc<-ClusDoc[!duplicated(ClusDoc)]


if(length(ClusPal)>0&length(ClusDoc)>0){
xl=c(min(Ccoor[ClusPal,axes[1]],Rcoor[ClusDoc,axes[1]]),max(Ccoor[ClusPal,axes[1]],Rcoor[ClusDoc,axes[1]]))
yl=c(min(Ccoor[ClusPal,axes[2]],Rcoor[ClusDoc,axes[2]]),max(Ccoor[ClusPal,axes[2]],Rcoor[ClusDoc,axes[2]]))
}

if(length(ClusPal)>0&length(ClusDoc)==0){
xl=c(min(Ccoor[ClusPal,axes[1]]),max(Ccoor[ClusPal,axes[1]]))
yl=c(min(Ccoor[ClusPal,axes[2]]),max(Ccoor[ClusPal,axes[2]]))
}

if(length(ClusPal)==0&length(ClusDoc)>0){
xl=c(min(Rcoor[ClusDoc,axes[1]]),max(Rcoor[ClusDoc,axes[1]]))
yl=c(min(Rcoor[ClusDoc,axes[2]]),max(Rcoor[ClusDoc,axes[2]]))
}

if(length(ClusPal)>0|length(ClusDoc)>0){
plot(0,0,pch=16,col="white",xlab=paste("DIM",axes[1],"(",round(eigen[axes[1],2],2),"%)",sep=""),ylab=paste("DIM",axes[2],"(",round(eigen[axes[2],2],2),"%)",sep=""),asp=1,main=main,xlim=xl,ylim=yl)
abline(h=0,v=0,lty=2)
if(length(ClusPal)>0){
points(Ccoor[ClusPal,axes[1]],Ccoor[ClusPal,axes[2]],col="red",pch=17)
text(Ccoor[ClusPal,axes[1]],Ccoor[ClusPal,axes[2]],rownames(Ccoor[ClusPal,]),col="red",cex=0.75,font=2,pos=3)
}
if(length(ClusDoc)>0){
points(Rcoor[ClusDoc,axes[1]],Rcoor[ClusDoc,axes[2]],col="blue",pch=16)
text(Rcoor[ClusDoc,axes[1]],Rcoor[ClusDoc,axes[2]],rownames(Rcoor[ClusDoc,]),col="blue",cex=0.75,font=2,pos=3)
}

}else{
print("There is no element to be represented")
}}
return (res=list( Metakeys.Metadocs=  Metakeys.MetaDocs, DimWord=DimensionWord))
}
