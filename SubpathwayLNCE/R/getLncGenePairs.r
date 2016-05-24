getLncGenePairs<-function(GeneExp,LncExp,a=0.025){
if(!exists("envData")) envData<-initialize()
pp<-get("pp",envir=envData)
code<-get("code",envir=envData)
nocode<-get("nocode",envir=envData)
GeneExp<-GeneExp[intersect(rownames(GeneExp),code),]
LncExp<-LncExp[intersect(rownames(LncExp),nocode),]
GeneExp<-ExpProcess(GeneExp,0.5)
LncExp<-ExpProcess(LncExp,0.5)
location<-matrix(0,length(pp[,1]),2)
location[,1]<-match(pp[,1],rownames(LncExp))
location[,2]<-match(pp[,2],rownames(GeneExp))
location<-na.omit(location)
dataset.1<-LncExp[location[,1],]
dataset.2<-GeneExp[location[,2],]
corValue<-c()
Pvalue<-c()
	   for(i in 1:length(location[,1])){
	       Pvalue[i]<-cor.test(as.numeric(dataset.1[i,]),as.numeric(dataset.2[i,]))$p.value
	       corValue[i]<-cor.test(as.numeric(dataset.1[i,]),as.numeric(dataset.2[i,]))$estimate
		  }   
LncGene<-cbind(rownames(dataset.1),rownames(dataset.2),corValue,Pvalue)
loca<-which(as.numeric(LncGene[,3])>0&as.numeric(LncGene[,4])<a)
LncGene.p<-LncGene[loca,c(1,2)]
colnames(LncGene.p)<-c("Lnc","Gene")
for(i in 1:length(LncGene.p[,1])){
 LncGene.p[i,1]<-unlist(strsplit(LncGene.p[i,1],".", fixed = T))[[1]]
 LncGene.p[i,2]<-unlist(strsplit(LncGene.p[i,2],".", fixed = T))[[1]]
}
mart<-get("mart",envir=envData)
lnc2Name<-get("lnc2Name",envir=envData)
LncGenePairs<-merge(LncGene.p,mart,by.x=colnames(LncGene.p)[2],by.y=colnames(mart)[1])[,c(2,3)]
LncGenePairs<-merge(LncGenePairs,lnc2Name,by.x=colnames(LncGenePairs)[1],by.y=colnames(lnc2Name)[1])[,c(3,2)]
return(LncGenePairs)

}
