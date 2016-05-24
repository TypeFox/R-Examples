Median.test<- function(y,trt,correct=TRUE,simulate.p.value = FALSE,console=TRUE){
name.y <- paste(deparse(substitute(y)))
name.t <- paste(deparse(substitute(trt)))
main<-paste(name.y,"~", name.t)
trt<-as.character(trt)
trt<-as.factor(trt)
junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
M0<-median(y)
medians <- tapply.stat(y,trt,stat="median") # change
above1<-tapply.stat(y,trt,function(x) sum(x>M0))
above2<-tapply.stat(y,trt,function(x) sum(x<=M0))
medians<-data.frame(medians,grather=above1[,2],lessEqual=above2[,2])
names(medians)[2]<-"Median"
ntr<-nrow(medians)
above<-y>M0
B<-table(above,trt)
O<-suppressWarnings(chisq.test(B,correct,simulate.p.value))
Tall<-O$statistic
Pall<-O$p.value
parameter<-O$parameter
# All comparison
comb <-utils::combn(ntr,2)
nn<-ncol(comb)
Median<-rep(0,nn)
Tstat<-rep(0,nn)
sig<-NULL
pvalue<-rep(0,nn)
for (k in 1:nn) {
i<-comb[1,k]
j<-comb[2,k]
ai<-medians[i,1]
aj<-medians[j,1]
subgroup<- junto[junto[,2]==ai | junto[,2]==aj,]
subgroup[,2]<-as.character(subgroup[,2])
subgroup[,2]<-as.factor(subgroup[,2])
M<-median(subgroup[,1])
above<-subgroup[,1]>M
B<-table(above,subgroup[,2])
O<-suppressWarnings(chisq.test(B,correct,simulate.p.value))
pvalue[k]<-O$p.value
Tstat[k]<-O$statistic
Median[k]<-M
sig[k]<-" "
if (pvalue[k] <= 0.001) sig[k]<-"***"
else  if (pvalue[k] <= 0.01) sig[k]<-"**"
else  if (pvalue[k] <= 0.05) sig[k]<-"*"
else  if (pvalue[k] <= 0.1) sig[k]<-"."
}
pvalue<-round(pvalue,4)
tr.i <- medians[comb[1, ],1]
tr.j <- medians[comb[2, ],1]
comparison<-data.frame(Median = Median,"Chisq"=Tstat, pvalue=pvalue,"sig"=sig)
rownames(comparison)<-paste(tr.i,tr.j,sep=" and ")
if(console){
cat("\nThe Median Test for",main,"\n")
cat("\nChi-square =", Tall,"  DF =", parameter,"  P.value", Pall)
cat("\nMedian =",M0,"\n\n")
print(comparison)
}
statistics<-data.frame(Chisq=Tall,p.chisq=Pall,Median=M0)
parameters<-data.frame(Df=parameter,ntr = ntr,test="Median")
rownames(parameters)<-" "
rownames(statistics)<-" "
output<-list(statistics=statistics,parameters=parameters,Medians=medians,
comparison=comparison,data=junto)
invisible(output)
}

