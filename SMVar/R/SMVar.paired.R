`SMVar.paired` <-
function(geneNumbers,logratio,fileexport=NULL,minrep=2,method="BH",threshold=0.05)
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
stopifnot(dim(geneNumbers)[1]==dim(logratio)[1])
nbgenes=dim(geneNumbers)[1]
globind=apply(logratio,1,FUN=function(x){sum(is.finite(x))>=minrep& var(x,na.rm=TRUE)>0})

globind=as.vector(which(globind==FALSE))

if (length(globind)>0)
{
logratio=logratio[-globind,]
geneNumbers=geneNumbers[-globind,]
nbgenes=nbgenes-length(globind)
print(paste(c("Warning:",(length(globind)),"gene(s) (is) are deleted because of too many missing values or null variance"),collapse=" "))
} 

nbrep=dim(logratio)[2]
ddl=nbrep-1
invddl=1/ddl

mean<-apply(logratio, 1, FUN = function(x) mean(x[is.finite(x)]))
RSS<-apply(logratio - mean, 1, FUN = function(x) sum(x[is.finite(x)]^2))

lnevari=log(RSS*invddl)
#datmixti=data.frame(geneId,lnevari)
#modeli=lm(lnevari~1,data=datmixti)
mui=1/length(lnevari)*sum(lnevari)
taui2=var(lnevari)-2/ddl

lambdai=taui2/(taui2+2/ddl)
sigma2=exp(mui+lambdai*(lnevari-mui))
varsigma2=sigma2^2*1/(1/taui2+ddl/2)

deltag=mean
teststat=deltag/sqrt(sigma2/nbrep)
ddlt=(2*sigma2^2)/varsigma2
Studentpval=2*(1-pt(abs(teststat),ddlt))
stat=data.frame(geneNumbers,teststat,Studentpval,ddlt,deltag)
colnames(stat)=c(colnames(geneNumbers),"TestStat","StudentPValue","DegOfFreedom","LogRatio")

n=length(colnames(geneNumbers))
adjpvalStud=p.adjust(stat$StudentPValue,method)
statStud=data.frame(stat,adjpvalStud)
colnames(statStud)=c(colnames(stat),"AdjPValue")
statStud=statStud[order(statStud$AdjPValue),]
genesdiffStud=statStud[which(statStud$AdjPValue<=threshold),c(1:length(colnames(geneNumbers)),(n+4),(n+5))]
print(paste(length(unique(genesdiffStud[,1])),"differentially expressed gene(s)",sep=" "))

if (!is.null(fileexport))
{
write.table(genesdiffStud[,-1],file=fileexport,sep="\t",row.names=FALSE)
}
else
{
print("Warning: No file is given to save the results of the analysis")
}
invisible(statStud)
}

