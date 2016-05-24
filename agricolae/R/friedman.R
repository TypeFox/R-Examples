`friedman` <-
function(judge,trt,evaluation,alpha=0.05,group=TRUE,main=NULL,console=FALSE){
name.x <- paste(deparse(substitute(judge)))
name.y <- paste(deparse(substitute(evaluation)))
name.t <- paste(deparse(substitute(trt)))
name.j <- paste(deparse(substitute(judge)))
if(is.null(main))main<-paste(name.y,"~", name.j,"+",name.t)
datos <- data.frame(judge, trt, evaluation)
matriz <- by(datos[,3], datos[,1:2], function(x) mean(x,na.rm=TRUE))
matriz <-as.data.frame(matriz[,])
#matriz <-as.matrix(evaluation)
name<-as.character(colnames(matriz))
ntr <-length(name)
m<-dim(matriz)
v<-array(0,m)
for (i in 1:m[1]){
v[i,]<-rank(matriz[i,])
}
vv<-as.numeric(v)
junto <- data.frame(evaluation, trt)
Means <- tapply.stat(junto[,1],junto[,2],stat="mean")  # change
sds <-   tapply.stat(junto[,1],junto[,2],stat="sd")    # change
nn <-   tapply.stat(junto[,1],junto[,2],stat="length") # change
mi<-tapply.stat(junto[,1],junto[,2],stat="min") # change
ma<-tapply.stat(junto[,1],junto[,2],stat="max") # change
nr<-unique(nn[,2])
s<-array(0,m[2])
# Suma de rangos por tratamiento
for (j in 1:m[2]){
s[j]<-sum(v[,j])
}
Means<-data.frame(Means,std=sds[,2],r=nn[,2],Min=mi[,2],Max=ma[,2]) 
names(Means)[1:2]<-c(name.t,name.y)
means<-Means[,c(1:2,4)]
rownames(Means)<-Means[,1]
Means<-Means[,-1]
means[,2]<-s
# row.names(means)<-means[,1]
rs<-array(0,m[2])
rs<-s-m[1]*(m[2]+1)/2
T1<-12*t(rs)%*%rs/(m[1]*m[2]*(m[2]+1))
T2<-(m[1]-1)*T1/(m[1]*(m[2]-1)-T1)
# Impresion de resultados
if(console){
cat("\nStudy:",main,"\n\n")
cat(paste(name.t,",",sep="")," Sum of the ranks\n\n")
print(data.frame(row.names = means[,1], means[,-1]))
cat("\nFriedman's Test")
cat("\n===============")
}
A1<-0
for (i in 1:m[1]) A1 <- A1 + t(v[i,])%*%v[i,]
DFerror <-(m[1]-1)*(m[2]-1)
Tprob<-qt(1-alpha/2,DFerror)
#
LSD<-as.numeric(Tprob*sqrt(2*(m[1]*A1-t(s)%*%s)/DFerror))
C1 <-m[1]*m[2]*(m[2]+1)^2/4
T1.aj <-(m[2]-1)*(t(s)%*%s-m[1]*C1)/(A1-C1)
T2.aj <-(m[1]-1)*T1.aj/(m[1]*(m[2]-1)-T1.aj)
p.value<-1-pchisq(T1.aj,m[2]-1)
p.noadj<-1-pchisq(T1,m[2]-1)
PF<-1-pf(T2.aj, ntr-1, (ntr-1)*(nr-1) )
if(console){
cat("\nAdjusted for ties")
cat("\nValue:",T1.aj)
cat("\nPvalue chisq :",p.value)
cat("\nF value :",T2.aj)
cat("\nPvalue F:",PF)
cat("\n\nAlpha     :",alpha)
cat("\nt-Student :",Tprob)
}
#...............
#cat("\nReplication:\t",nr)
if (group) {
if(console){
cat("\nLSD       :",LSD)
cat("\n\nMeans with the same letter are not significantly different.")
cat("\nGroupTreatment and Sum of the ranks\n")}
s<-as.numeric(s)
groups<-order.stat(name,s,LSD,console=console)
names(groups)[2]<-"Sum of ranks"
comparison=NULL
statistics<-data.frame(Chisq=T1.aj,p.chisq=p.value,F=T2.aj,p.F=PF,LSD)
}
if (!group) {
comb <-utils::combn(ntr,2)
nn<-ncol(comb)
dif<-rep(0,nn)
pvalue<-rep(0,nn)
LCL<-dif
UCL<-dif
sig<-NULL
LSD<-rep(0,nn)
stat<-rep("ns",nn)
for (k in 1:nn) {
i<-comb[1,k]
j<-comb[2,k]
#if (means[i, 2] < means[j, 2]){
#comb[1, k]<-j
#comb[2, k]<-i
#}
dif[k]<-s[comb[1,k]]-s[comb[2,k]]
sdtdif<- sqrt(2*(m[1]*A1-t(s)%*%s)/DFerror)
pvalue[k]<- round(2*(1-pt(abs(dif[k])/sdtdif,DFerror)),4)
LSD[k]<-round(Tprob*sdtdif,2)
LCL[k] <- dif[k] - LSD[k]
UCL[k] <- dif[k] + LSD[k]
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
if(console){cat("\n\nComparison between treatments\nSum of the ranks\n\n")
print(comparison)}
statistics<-data.frame(Chisq=T1.aj,p.chisq=p.value,F=T2.aj,p.F=PF)
groups=NULL
# output<-data.frame(trt= means[,1],means= means[,2],M="",N=means[,3])
}
parameters<-data.frame(Df=ntr-1,ntr = ntr, t.value=Tprob,alpha=alpha,test="Friedman",name.t=name.t)
rownames(parameters)<-" "
rownames(statistics)<-" "
Means<-data.frame(rankSum=means[,2],Means)
output<-list(statistics=statistics,parameters=parameters, 
		means=Means,comparison=comparison,groups=groups)
invisible(output)
}
