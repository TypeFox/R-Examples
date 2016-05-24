`durbin.test` <-
function(judge,trt,evaluation,alpha=0.05, group=TRUE,main=NULL,console=FALSE) {
name.y <- paste(deparse(substitute(evaluation)))
name.t <- paste(deparse(substitute(trt)))
name.j <- paste(deparse(substitute(judge)))
if(is.null(main))main<-paste(name.y,"~", name.j,"+",name.t)
judge<-as.factor(judge)
trt<-as.factor(trt)
k <-unique(table(judge))
r <-unique(table(trt))
b <-nlevels(judge)
ntr <-nlevels(trt)
lambda<-r*(k-1)/(ntr-1)
x<-data.frame(judge,trt,evaluation)
Means<- tapply.stat(x[,3],x[,2],stat="mean")  # change
sds <-  tapply.stat(x[,3],x[,2],stat="sd")
mi <-   tapply.stat(x[,3],x[,2],stat="min")
ma <-   tapply.stat(x[,3],x[,2],stat="max")
Means<-data.frame(Means,std=sds[,2],r,Min=mi[,2],Max=ma[,2])
rownames(Means)<-Means[,1]
Means<-Means[,-1]
names(Means)[1] <- name.y
# Determina el rango dentro de cada juez
z <- by(x,x$judge,function(x) rank(x$evaluation))
y<-data.frame(c(z))
m<-dim(y)
n<-m[1]*m[2]
rango <- 1:n
for (i in 1:m[1]) {
for (j in 1:m[2]) {
kk=i+m[1]*(j-1)
rango[kk]<-y[i,j]
}
}
x <- data.frame(x, rango)
means <- tapply.stat(x[, 4], x[, 2], stat = "sum")
names(means)[1:2] <- c(name.t, name.y)
x<-data.frame(x,rango)

names(means)[1:2]<-c(name.t,name.y)
z <-by(x,x$trt,function(x) sum(x$rango))
y<-as.vector(c(z))
name<-as.character(dimnames(z)$"x$trt")
s <- (y-r*(k+1)/2)^2
s1 <- sum(s)
# determina el valor de Durbin
gl1<-ntr-1 ;gl2<-b*k-ntr-b+1
C <- b*k*(k+1)^2/4
A <- sum(rango^2)
s <- (ntr - 1) * s1/(A-C)
prob<-1-pchisq(s,gl1); Tprob<-qt(1-alpha/2,gl2)
sdtdif <- sqrt(2*r*(A-C)*(1-s/(b*(k-1)))/gl2)
LSD <-Tprob*sdtdif
nameTrt<-as.character(means[,1])
# s,prob,Tprob,Mc,gl1,gl2)
# Impresion de resultados
if(console){
cat("\nStudy:",main,"\n")
cat(paste(name.t,",",sep="")," Sum of ranks\n\n")
print(data.frame(row.names = nameTrt, sum=means[,2]))
cat("\nDurbin Test")
cat("\n===========")
cat("\nValue      :",s)
cat("\nDf 1       :",gl1)
cat("\nP-value    :",prob)
cat("\nAlpha      :",alpha)
cat("\nDf 2       :",gl2)
cat("\nt-Student  :",Tprob)
cat("\n\nLeast Significant Difference\nbetween the sum of ranks: ",LSD,"\n")
# comparacion de tratamientos.
cat("\nParameters BIB")
cat("\nLambda     :",lambda)
cat("\ntreatmeans :",ntr)
cat("\nBlock size :",k)
cat("\nBlocks     :",b)
cat("\nReplication:",r,"\n")
}
if (group)
{
if(console)cat("\nGroups, Treatments and sum of the ranks\n\n")
y<-as.numeric(y)
groups<-order.stat(name,y,LSD,console=console)
comparison<-NULL
}

if (!group) {
comb <-utils::combn(ntr,2)
nn<-ncol(comb)
dif<-rep(0,nn)
pvalue<-rep(0,nn)
sig<-rep(" ",nn)
for (kk in 1:nn) {
i<-comb[1,kk]
j<-comb[2,kk]
#if (y[i] < y[j]){
#comb[1,kk]<-j
#comb[2,kk]<-i
#}
dif[kk]<-y[comb[1,kk]]-y[comb[2,kk]]
pvalue[kk]<- 2*round(1-pt(abs(dif[kk])/sdtdif,gl2),4)
sig[kk]<-" "
if (pvalue[kk] <= 0.001) sig[kk]<-"***"
else  if (pvalue[kk] <= 0.01) sig[kk]<-"**"
else  if (pvalue[kk] <= 0.05) sig[kk]<-"*"
else  if (pvalue[kk] <= 0.1) sig[kk]<-"."
}
tr.i <- nameTrt[comb[1, ]]
tr.j <- nameTrt[comb[2, ]]
comparison<-data.frame("Difference" = dif, pvalue=pvalue,"sig."=sig)
rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
if(console){cat("\nComparison between treatments sum of the ranks\n\n")
print(comparison)}
groups=NULL
}
#output<-data.frame(means,M="",N=r)
#
parameters<-data.frame(lambda=lambda,treatments=ntr,blockSize=k,blocks=b,r=r,alpha=alpha,test="Durbin",name.t=name.t)
statistics<-data.frame(chisq.value=s, p.value=prob, t.value=Tprob,LSD=LSD)
	rownames(parameters)<-" "
	rownames(statistics)<-" "
	output<-list(statistics=statistics,parameters=parameters, 
	means=Means,rank=means,comparison=comparison,groups=groups)
	
invisible(output)
}
