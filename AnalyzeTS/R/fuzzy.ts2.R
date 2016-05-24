fuzzy.ts2 <-
function(ts,n=5,w=NULL,D1=0,D2=0,C=NULL,r=4,trace=FALSE,forecast=NULL,plot=FALSE){

#---ham con---
is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  
abs(x - round(x)) < tol

namthang<-function(data.ts){
batdau<-start(data.ts)
tanso<-frequency(data.ts)
nam1<-batdau[1]
thang1<-batdau[2]
ketthuc<-end(data.ts)
nam2<-ketthuc[1]
thang2<-ketthuc[2]
namkq<-1:length(data.ts)
thangkq<-1:length(data.ts)
index=0;
for(nam in nam1:nam2)
for(thang in 1:tanso)
if(nam!=nam1 || thang>=thang1){
index=index+1
namkq[index]<-nam
thangkq[index]<-thang
if(nam==nam2 & thang==thang2) break;
}
if(tanso==4)
{thangkq[thangkq==1]<-"Q1";thangkq[thangkq==2]<-"Q1";
thangkq[thangkq==3]<-"Q1";thangkq[thangkq==4]<-"Q1";
print<-paste(namkq,thangkq,sep=" ")}
else
if(tanso==12)
{thangkq[thangkq==1]<-"Jan";thangkq[thangkq==2]<-"Feb";
thangkq[thangkq==3]<-"Mar";thangkq[thangkq==4]<-"Apr";
thangkq[thangkq==5]<-"May";thangkq[thangkq==6]<-"Jun";
thangkq[thangkq==7]<-"Jul";thangkq[thangkq==8]<-"Aug";
thangkq[thangkq==9]<-"Sep";thangkq[thangkq==10]<-"Oct";
thangkq[thangkq==11]<-"Nov";thangkq[thangkq==12]<-"Dec";
print<-paste(namkq,thangkq,sep=" ")}
else
if(tanso==7)
{thangkq[thangkq==1]<-"Mon";thangkq[thangkq==2]<-"Tue";
thangkq[thangkq==3]<-"Wed";thangkq[thangkq==4]<-"Thu";
thangkq[thangkq==5]<-"Fri";thangkq[thangkq==6]<-"Sat";
thangkq[thangkq==7]<-"Sun";
print<-paste(namkq,thangkq,sep=" ")}
else
if(tanso!=1)
print<-paste("(",namkq,",",thangkq,")",sep="")
else
print<-namkq
print
}

computeVi<-function(M,table,w){
n<-(dim(M)[1] - w)
cot<-dim(M)[2]
Vi<-1:dim(M)[1]
Vi[1:w]<-NA
for(i in 1:n){
O<-M[i:(w-1+(i-1)),]
if(w==2)O<-t(as.matrix(O))
K<-M[w+(i-1),]
R<-O

for (i1 in 1:cot) for (j1 in 1:(w-1)) if (O[j1,i1] > K[i1]) R[j1,i1] <- K[i1]
F <- 1:cot
for (i2 in 1:cot) F[i2] <- max(R[, i2])
Vi[w+i]<- sum(F * table$Bw)/sum(F)
}
Vi
}

computeVi2<-function(M,table){
cot<-dim(M)[2]
dong<-(dim(M)[1]-1)
O<-M[1:dong,]
if(w==2)O<-t(as.matrix(O))
K<-M[(dong+1),]
R<-O
for (i1 in 1:cot) for (j1 in 1:(w-1)) if (O[j1,i1] > K[i1]) R[j1,i1] <- K[i1]
F <- 1:cot
for (i2 in 1:cot) F[i2] <- max(R[, i2])
Vi<- sum(F * table$Bw)/sum(F)

Vi
}
#---ham con---


#---kiem tra input---
if (!is.numeric(ts)) stop("Error in 'ts'!")
if(!is.ts(ts)) stop("Error in 'ts'!") else if(!is.null(dim(ts)))stop("Error in 'ts'!")

kt<-0
for(i in 1:length(ts))if(is.na(ts[i]))kt=kt+1
if(kt>0)stop("Time series contain NA!")

if (is.na(n)||!is.numeric(n) || n < 1 || !is.wholenumber(n)) stop("Error in 'n'!")
if (is.na(r)||!is.numeric(r) || r < 0 || !is.wholenumber(r)) stop("Error in 'r'!")
if (is.null(w)||is.na(w)||!is.numeric(w) || w < 2 || !is.wholenumber(w)) stop("Error in 'w'!")
if (is.null(C)||is.na(C)||!is.numeric(C)) stop("Error in 'C'!")
if (is.na(D1)||!is.numeric(D1)) stop("Error in 'D1'!")
if (is.na(D2)||!is.numeric(D2)) stop("Error in 'D2'!")
if (is.null(forecast)||is.na(forecast)||!is.numeric(forecast) || forecast < 1 || !is.wholenumber(forecast)) stop("Error in 'forecast'!")
if(w>=length(ts)) stop("Error in 'w'!")
#---kiem tra input---


#---tap bien doi---
ts1<-as.vector(diff(ts))
min.x=min(ts1)-D1
max.x=max(ts1)+D2
h=(max.x-min.x)/n
k<-1:(n+1);U<-1:n

for(i in 1:(n+1)){
if(i==1)k[i]=min.x
else
{k[i]=min.x+(i-1)*h
U[i-1]=paste("u",i-1,sep="")}}

D<-data.frame(U,low=k[1:n],up=k[2:(n+1)])
D$Bw<-(1/2)*(D$low+D$up)
table1<-D
#---tap bien doi---


#---mo hoa bien doi---
thoidiem <- namthang(ts)
Ai <- 1:length(ts1)
MATRIX <- matrix(1:(length(ts1) * n), ncol = n)
for (i in 1:length(ts1)) {
temp = ""
for (j in 1:n) {
my.At <- 1/(1 + (C * (ts1[i] - D$Bw[j]))^2)
MATRIX[i, j] <- my.At
At.j <- paste("(",round(my.At,r),"/u",j,sep = "", ")")
if (j == 1) temp <- paste(temp, At.j, sep = "")
else temp <- paste(temp, At.j, sep = ",")
}
Ai[i] <- paste("A[",thoidiem[i+1],"]={",temp,"}",sep="")
}
table2<-data.frame(point=thoidiem,ts=ts,diff.ts=c(NA,ts1))
table3<-c(NA,Ai)
#---mo hoa bien doi---


#---noi suy---
V<-computeVi(MATRIX,table1,w)
N<-c(NA,(ts[-length(ts)]+V))
V<-c(NA,V)
table4<-data.frame(point=thoidiem,interpolate=N,diff.interpolate=V)
accuracy<-av.res(Y=data.frame(ts),F=data.frame("Abbasov.Mamedova"=table4[,2]))
table4<-na.omit(table4)
rownames(table4)<-c(1:dim(table4)[1])
#---noi suy---


#---du bao---
V<-1:forecast
N<-0:forecast
Ai<-1:forecast
N[1]<-table4[dim(table4)[1],2]

MT<-MATRIX[(dim(MATRIX)[1]-w+1):dim(MATRIX)[1],]

temp<-c(as.vector(ts),1:forecast)
temp<-ts(temp,start=start(ts),frequency=frequency(ts))
temp<-namthang(temp)
temp<-temp[(length(ts)+1):length(temp)]
thoidiem<-temp

for(i in 1:forecast){
V[i]<-computeVi2(MT,table1)
N[i+1]<-N[i]+V[i]

for(chuyen in 1:(dim(MT)[1]-1)) MT[chuyen,]<-MT[(chuyen+1),]

temp = ""
for (j in 1:n) {
my.At <- 1/(1 + (C * (V[i] - D$Bw[j]))^2)
MT[dim(MT)[1], j] <- my.At

At.j <- paste("(",round(my.At,r),"/u",j,sep = "", ")")
if (j == 1) temp <- paste(temp, At.j, sep = "")
else temp <- paste(temp, At.j, sep = ",")
}
Ai[i] <- paste("A[",thoidiem[i],"]={",temp,"}",sep="")
}

N<-N[-1]
table5<-data.frame(point=thoidiem,forecast=N,diff.forecast=V)
table6<-Ai
#---du bao---

KQ<-list(type="Abbasov-Manedova",table1=table1,table2=table2,table3=table3,table4=table4,table5=table5,table6=table6,accuracy=accuracy)
if(trace==TRUE) MO<-KQ else if(trace==FALSE) MO<-table5 else MO<-c("'trace' must be 'TRUE' or 'FALSE'")

if(plot==TRUE) {
tsp<-ts
tsp[1:w]<-NA
tsp<-na.omit(tsp)

goc<-ts(table2[,2],start=start(ts),frequency=frequency(ts))
dubao<-ts(c(table4[,2],table5[,2]),start=start(tsp),frequency=frequency(tsp))

layout(1:1)
ts.plot(goc)

if(length(goc)<50){
plot(goc,col="blue",main=paste("Abbasov-Mamedova, C =",C,"w =",w,"n =",n),type="o",
     pch=16,ylim=c(min(c(goc,dubao)),max(c(goc,dubao))),xlim=c(par()$xaxp[1],(par()$xaxp[2]+(forecast+2)/frequency(ts))),
     xlab="index",ylab="data");
lines(dubao,col="red",type="o",pch=18)
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1),pch=c(16,18))
}

if(length(goc)>49){
plot(goc,col="blue",main=paste("Abbasov-Mamedova, C =",C,"w =",w,"n =",n),type="l",
     ylim=c(min(c(goc,dubao)),max(c(goc,dubao))),xlim=c(par()$xaxp[1],(par()$xaxp[2]+(forecast+5)/frequency(ts))),
     xlab="index",ylab="data");
lines(dubao,col="red",type="l")
legend("bottomright","(x,y)",c("ts","forecast"),col=c("blue","red"),lty=c(1,1))
}

k.ve<-par()$yaxp
h=(k.ve[2]-k.ve[1])/k.ve[3]
for(i in -2:(k.ve[3]+2)) abline(h=k.ve[1]+h*i,lty=3,col="gray")

}


MO
}
