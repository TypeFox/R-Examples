`waerden.test` <-
function(y, trt,alpha=0.05,group=TRUE,main=NULL,console=FALSE) {
name.y <- paste(deparse(substitute(y)))
name.t <- paste(deparse(substitute(trt)))
if(is.null(main))main<-paste(name.y,"~", name.t)
junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
Means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
sds <-   tapply.stat(junto[,1],junto[,2], stat="sd")  #change
nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
mi<-tapply.stat(junto[,1],junto[,2],stat="min") # change
ma<-tapply.stat(junto[,1],junto[,2],stat="max") # change
Means<-data.frame(Means,std=sds[,2],r=nn[,2],Min=mi[,2],Max=ma[,2]) 
rownames(Means)<-Means[,1]
Means<-Means[,-1]
names(Means)[1]<-name.y

N<- nrow(junto)
junto[, 1] <- qnorm(round(rank(junto[, 1]) /(N+1),3))
S <- sum(junto[,1]^2)/(N-1)
means <- tapply.stat(junto[,1],junto[,2],stat="mean") # change
sds <-   tapply.stat(junto[,1],junto[,2], stat="sd")  #change
nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
means<-data.frame(means,std=sds[,2],r=nn[,2])  
names(means)[1:2]<-c(name.t,name.y)
#row.names(means)<-means[,1]
ntr<-nrow(means)
DFerror<-N - ntr
T1 <- 0
for (i in 1:ntr) {
T1 <- T1 + means[i, 2]^2*means[i,4] # change
}
T1<-T1/S
p.chisq <- 1 - pchisq(T1, ntr - 1)
if(console){
cat("\nStudy:",main)
cat("\nVan der Waerden (Normal Scores) test's\n")
cat("\nValue :", T1)
cat("\nPvalue:", p.chisq)
cat("\nDegrees of freedom: ", ntr - 1,"\n\n")
cat(paste(name.t,",",sep="")," means of the normal score\n\n")
print(data.frame(row.names = means[,1], means[,-1]))
}
MSerror <- S * ((N - 1 - T1)/(N - ntr))
#cat("\nComparison of treatments")
#...............

nr <- unique(means[,4]) # change
nr1<-nr
Tprob<-qt(1-alpha/2,DFerror)
if (group) {
if(console){cat("\nt-Student:", Tprob)
cat("\nAlpha    :",alpha)}
    if (length(nr1) == 1) {
        LSD <- Tprob * sqrt(2 * MSerror/nr)
if(console)cat("\nLSD      :", LSD,"\n")
statistics<-data.frame(Chisq=T1,p.chisq=p.chisq,LSD=LSD )
    }
    else {
    if(console)cat("\nMinimum difference changes for each comparison\n")
#         nr <- 1/mean(1/nn[, 2])
#         LSD <- Tprob * sqrt(2 * MSerror/nr)
#         cat("\nLSD      :", LSD,"\n")
#         cat("\nHarmonic Mean of Cell Sizes ", nr)
		 statistics<-data.frame(Chisq=T1,p.chisq=p.chisq)
	 }   
if(console){cat("\nMeans with the same letter are not significantly different\n")
cat("\nGroups, Treatments and means of the normal score\n")}
groups <- order.group(means[,1], means[,2], means[,4], MSerror, Tprob,std.err=sqrt(MSerror/ means[,4]),console=console) # change
groups<-groups[,1:3]
comparison=NULL
 }
 if (!group) {
comb <-utils::combn(ntr,2)
nn<-ncol(comb)
dif<-rep(0,nn)
LCL<-dif
UCL<-dif
sig<-NULL
pvalue<-rep(0,nn)
for (k in 1:nn) {
i<-comb[1,k]
j<-comb[2,k]
#if (means[i, 2] < means[j, 2]){
#comb[1, k]<-j
#comb[2, k]<-i
#}
dif[k]<-means[i,2]-means[j,2]
sdtdif<- sqrt(S*((N-1-T1)/(N-ntr))*(1/means[i,4]+1/means[j,4])) # change
pvalue[k]<- round(2*(1-pt(abs(dif[k])/sdtdif,DFerror)),4)
LSD <- Tprob*sdtdif
LCL[k] <- dif[k] - LSD
UCL[k] <- dif[k] + LSD
sig[k]<-" "
if (pvalue[k] <= 0.001) sig[k]<-"***"
else  if (pvalue[k] <= 0.01) sig[k]<-"**"
else  if (pvalue[k] <= 0.05) sig[k]<-"*"
else  if (pvalue[k] <= 0.1) sig[k]<-"."
}
tr.i <- means[comb[1, ],1]
tr.j <- means[comb[2, ],1]
comparison<-data.frame("Difference" = dif, pvalue=pvalue,"sig."=sig,LCL,UCL)
rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
if(console){cat("\nComparison between treatments means\nmean of the normal score\n\n")
print(comparison)}
statistics<-data.frame(Chisq=T1,p.chisq=p.chisq) 
groups=NULL
#output<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,4]) # change
}
Means<-data.frame(normlScore=means[,2],Means)
parameters<-data.frame(Df=ntr-1,ntr = ntr, t.value=Tprob,alpha=alpha,test="Waerden",name.t=name.t)
rownames(parameters)<-" "
rownames(statistics)<-" "
output<-list(statistics=statistics,parameters=parameters, 
		means=Means,comparison=comparison,groups=groups)
    invisible(output)
}

