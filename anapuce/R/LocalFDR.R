LocalFDR <-
function(dataf = dataf, graph = TRUE, method = NULL,lambda0 = 0.5, smoothing="1", thres=c(0.01,0.05,0.10,0.20), mm=c(3,5,15,NA)) {

# We ordered the pvalues
dataord<-dataf[order(as.numeric(dataf[,2])),]
p<-dataord[,2]
n<-length(p)-2

# Calculation of m0 with lambda0=0.5 by default

listemethode<-c("conservative","adaptive","bootstrap","smoother")

if ((method %in% listemethode)==TRUE && (is.null(method) == FALSE)) {
p<-as.numeric(p)
# We used the function fdr.estimate.eta0 from the package GeneTS
#(written by Konstantinos Fokianos and Korbinian Strimmer and   
# adapted in part from S-PLUS code by Y. Benjamini and R code from J.D. Storey. 

m0 <- fdr.estimate.eta0(p,method=method)*length(p)
} else {
m0 <- 2*sum((p>=lambda0)*1)
p<-as.numeric(p)
}


if (smoothing=="1"){
# Step 1: Calculation of FDRU and FDRV

tmpu <- diff(p)
u<-tmpu[1:(length(tmpu)-1)]
v<-diff(p,2)/2

FDRu<-cbind(2:(n+1),m0*u)
FDRv<-cbind(3:(n+2),m0*v)

# Step 2: Smooth the FDRu and FDRv curves in function of i using a lowess function (f=0.2)

FDRufit<-loess(FDRu[,2]~FDRu[,1],span=0.2)$fitted
FDRvfit<-loess(FDRv[,2]~FDRv[,1],span=0.2)$fitted

# Moving average smoothing

# Remark : The number of intervals goes from 1 to 5
# Convention :  If mm[4]=NA then mm[4]=m0/4
if (is.na(mm[4])) mm[4]<-round(m0/4)

indice_u <- findInterval(FDRu[,2],thres)+1
FDRuL<-apply(as.matrix(FDRu[,2]),2,lissageMM,indice=indice_u,mm)
indice_v <- findInterval(FDRv[,2],thres)+1
FDRvL<-apply(as.matrix(FDRv[,2]),2,lissageMM,indice=indice_v,mm)

if (graph){
pdf(file=paste("LocalFDRGraph.pdf",sep=""))
plot(FDRv[,1],FDRv[,2],main="",xlab="",ylab="",type="l")
plot(FDRv[,1],FDRvL,main="",xlab="",ylab="",type="l",ylim=c(0,1.2))
lines(FDRv[,1],FDRvfit,type="l",col=2)
dev.off()
}
FDRuL <- pmin(c(0,FDRuL,1),1)
FDRvL <- pmin(c(0,FDRvL,1),1)
out<-cbind(dataord,c("NA",FDRu[,2],"NA"),c("NA",FDRv[,2],"NA"),c("NA",FDRufit,"NA"),c("NA",FDRvfit,"NA"),FDRuL,FDRvL)
colnames(out)<-c("gene","pvalue","FDRu","FDRv","smoothed_lowess_FDRu","smoothed_lowess_FDRv","smoothed_MA_FDRu","smoothed_MA_FDRv")
# Step 3 : Output data file


} else {

if (smoothing=="2") {
    diff.p <- diff(p)
    FDRuL <- m0 * isoreg(c(p[1], 
        diff.p))$yf
    FDRuL <- pmin(FDRuL, 1)
    diff.pv <- diff(p,2)/2
   # FDRvL <- pmin(c(0,m0 * isoreg(c(p[1], 
   #     diff.pv))$yf),1)
        
# Step 3 : Output data file
out<-cbind(dataord,FDRuL)
colnames(out)<-c("gene","pvalue","smoothed_PAVA_FDR")
    
}
}

write.table(out,file="LocalFDRFile.txt",quote=TRUE,row.names=FALSE,sep="\t")
invisible(out)
#  (c) 2007 Institut National de la Recherche Agronomique

}

