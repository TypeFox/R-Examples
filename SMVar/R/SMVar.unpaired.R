`SMVar.unpaired` <-
function(geneNumbers,listcond,fileexport=NULL,minrep=2,method="BH",threshold=0.05)
{
if (minrep<2){print("warning: minrep must be >=2")}
stopifnot(class(geneNumbers) %in% c("vector","matrix","data.frame","character","numeric","integer"))
if(class(geneNumbers) %in% c("matrix","data.frame"))
{
if (is.null(colnames(geneNumbers))){print("Error: geneNumbers must have column names")}
geneNumbers=data.frame(GeneId=rep(1:dim(geneNumbers)[1],1),geneNumbers)
}
if(class(geneNumbers) %in% c("vector","character","numeric","integer"))
{
geneNumbers=data.frame(GeneId=rep(1:length(geneNumbers),1),geneNumbers)
}
nbcond=length(listcond)
stopifnot(dim(geneNumbers)[1]==dim(listcond[[1]])[1])
nbgenes=dim(geneNumbers)[1]
for (i in 1:nbcond)
{
condi=listcond[[i]]
indcond=apply(condi,1,FUN=function(x){sum(is.finite(x))>=minrep& var(x,na.rm=TRUE)>0})
assign(paste("indcond",i,sep=""),indcond)
}

globind=get(paste("indcond",1,sep=""))
for (i in 2:nbcond)
{globind=globind & get(paste("indcond",i,sep=""))}

globind=as.vector(which(globind==FALSE))

if (length(globind)>0)
{
for (i in 1:nbcond)
{
listcond[[i]]=listcond[[i]][-globind,]
}
geneNumbers=geneNumbers[-globind,]
nbgenes=nbgenes-length(globind)
print(paste(c("Warning:",(length(globind)),"gene(s) (is) are deleted because of too many missing values or null variance"),collapse=" "))
} 

for (i in 1:nbcond)
{
condi=listcond[[i]]
nbreptemp=apply(condi,1,FUN=function(x) sum(is.finite(x)))
assign(paste("nbrep",i,sep=""),nbreptemp)
assign(paste("ddl",i,sep=""),(nbreptemp-1))
assign(paste("invddl",i,sep=""),1/get(paste("ddl",i,sep="")))
}

calcmRSS=function(datacond)
{
mean<-apply(datacond, 1, FUN = function(x) mean(x[is.finite(x)]))
RSS<-apply(datacond - mean, 1, FUN = function(x) sum(x[is.finite(x)]^2))
datacond<-cbind(mean,RSS,datacond)
datacond
}
listcondRSS=lapply(listcond,calcmRSS)
for (i in 1:nbcond)
{
lnevari=log(listcondRSS[[i]][,2]*get(paste("invddl",i,sep="")))
#datmixti=data.frame(geneId,lnevari)
#modeli=lm(lnevari~1,data=datmixti)
mui=1/length(lnevari)*sum(lnevari)
taui2=var(lnevari)-2/get(paste("ddl",i,sep=""))
lambdai=taui2/(taui2+2/get(paste("ddl",i,sep="")))
assign(paste(c("sigma2",i),collapse=""),exp(mui+lambdai*(lnevari-mui)))
assign(paste(c("varsigma2",i),collapse=""),(get(paste(c("sigma2",i),
collapse="")))^2*1/(1/taui2+get(paste("ddl",i,sep=""))/2))
}
for(j in 1:(nbcond-1))
{
for(k in (j+1):nbcond)
{
assign(paste(c("deltag",j,k),collapse=""),(listcondRSS[[k]][,1]-listcondRSS[[j]][,1]))
s2jtemp=get(paste(c("sigma2",j),collapse=""))
s2ktemp=get(paste(c("sigma2",k),collapse=""))
njtemp=get(paste("nbrep",j,sep=""))
nktemp=get(paste("nbrep",k,sep=""))
teststat=get(paste(c("deltag",j,k),collapse=""))/sqrt(s2jtemp/njtemp+s2ktemp/nktemp)
ddlt=(s2jtemp/njtemp+s2ktemp/nktemp)^2/
(get(paste(c("varsigma2",j),collapse=""))/(2*njtemp^2)+get(paste(c("varsigma2",k),collapse=""))/(2*nktemp^2))
Studentpval=2*(1-pt(abs(teststat),ddlt))
cond1=rep(j,nbgenes)
cond2=rep(k,nbgenes)
stat=data.frame(geneNumbers,teststat,Studentpval,ddlt,get(paste(c("deltag",j,k),collapse="")),cond1,cond2)
colnames(stat)=c(colnames(geneNumbers),"TestStat","StudentPValue","DegOfFreedom","LogRatio","Cond1","Cond2")
assign(paste(c("stat",j,k),collapse=""),stat)
}
}
globstat=data.frame()
for(j in 1:(nbcond-1))
{
for(k in (j+1):nbcond)
{
globstat=rbind(globstat,get(paste(c("stat",j,k),collapse="")))
}
}
n=length(colnames(geneNumbers))
adjpvalStud=p.adjust(globstat$StudentPValue,method)
globstatStud=data.frame(globstat,adjpvalStud)
colnames(globstatStud)=c(colnames(globstat),"AdjPValue")
globstatStud=globstatStud[order(globstatStud$AdjPValue),]
genesdiffStud=globstatStud[which(globstatStud$AdjPValue<=threshold),c(1:length(colnames(geneNumbers)),(n+4):(n+7))]
print(paste(length(unique(genesdiffStud[,1])),"differentially expressed gene(s)",sep=" "))

if (!is.null(fileexport))
{
write.table(genesdiffStud[,-1],file=fileexport,sep="\t",row.names=FALSE)
}
else
{
print("Warning: No file is given to save the results of the analysis")
}
invisible(globstatStud)
}

