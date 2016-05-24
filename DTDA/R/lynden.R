lynden <-
function(X, U=NA, V=NA, error=NA, nmaxit=NA,
 boot=TRUE, B=NA, alpha=NA, display.F=FALSE, 
display.S=FALSE){



trunc<-"double"

if(all(is.na(U))==TRUE){

trunc<-"right"
cat("case U=NA","\n")
if(any(is.na(V))==TRUE|any(is.na(X))==TRUE){
navec<-c(which(is.na(X)),which(is.na(V)))
X<-X[-navec]
V<-V[-navec]
}}
if(all(is.na(V))==TRUE){

trunc<-"left"
cat("case V=NA","\n")
if(any(is.na(U))==TRUE|any(is.na(X))==TRUE){
navec<-c(which(is.na(X)),which(is.na(U)))
X<-X[-navec]
U<-U[-navec]
}}

if(trunc=="double"){
if(any(is.na(U))==TRUE | any(is.na(V))==TRUE|any(is.na(X))==TRUE){

navec<-c(which(is.na(X)),which(is.na(U)),which(is.na(V)))
X<-X[-navec]
U<-U[-navec]
V<-V[-navec]

}}

D<-cbind(X,U,V)

if(all(is.na(U))==TRUE){
D[,2]<--D[,3]
D[,1]<--D[,1]
D[,3]<-rep(max(D[,1])+1,length(X))
}
if(all(is.na(V))==TRUE){
D[,3]<-rep(max(D[,1])+1,length(X))
}

ord<-order(D[,1])
C<-matrix(0,nrow=nrow(D),ncol=ncol(D))
C[,1]<-sort(D[,1])
C[,2:ncol(D)]<-D[ord,2:ncol(D)]


NJ<-matrix(data=0,ncol=1,nrow=nrow(C))
K<-NJ
for(j in 1:nrow(C)){
for(k in 1:nrow(C)){
if (C[k,2]<=C[j,1]& C[k,1]>=C[j,1]) {NJ[j,1]<-NJ[j,1]+1}
}}
K[,1]<-rep(1,length=nrow(C))


mNJ<-min(NJ[1:nrow(C)-1,1])
if (mNJ!=1) K[,1]<-NJ[,1]
if (mNJ==1){
istar<-max(which(NJ[1:nrow(C)-1,1]==1))
K[1:istar,1]<-2
if(istar<nrow(C)){K[(istar+1):nrow(C),1]<-NJ[(istar+1):nrow(C),1]}
}
NJ<-as.matrix(apply(cbind(NJ,K),1,max),ncol=1)


J<-matrix(data=0,ncol=nrow(C),nrow=nrow(C))
for (k in 1:nrow(C)){
for(j in 1:nrow(C)) {
a1<-min(C[j,2],C[j,3])
b1<-max(C[j,2],C[j,3])
if(C[k,1]>=a1 & C[k,1]<=b1) J[k,j]<-1
}}
h0<-1/NJ


G0<-matrix(data=0,ncol=1,nrow=nrow(C))
G0[1,1]<-1
for(j in 2:nrow(C)){
G0[j,1]<-exp(sum(log(1-h0[1:(j-1),1])))
}


f0<-matrix(data=G0[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f0[j,1]<-(G0[j,1]-G0[j+1,1])
}

Gvi0<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxi<-0
for (k in 1:nrow(C)){
if (C[j,3]< C[k,1]) {auxi<-sum(f0[k,])+ auxi}
}
Gvi0[j,]<-auxi
}

F0<-t(J)%*%f0
Q0<-Gvi0/F0
h<-h0
for(i in 1:nrow(h0)){
h[i,]<- 1/(NJ[i,]+(J[i,])%*%Q0)
}

S0<-1
if(is.na(error)==TRUE) error<-1e-6

if(is.na(nmaxit)==TRUE)nmaxit<-100
iter<-0

while(S0>error|iter>nmaxit){
iter<-iter+1
if (iter>nmaxit) stop("Default number of iterations not enough for convergence")
G1<-matrix(data=0,ncol=1,nrow=nrow(C))
G1[1,1]<-1
for(j in 2:nrow(C)){
G1[j,1]<-exp(sum(log(1-h[1:(j-1),1])))
}

f<-matrix(data=G1[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f[j,1]<-(G1[j,1]-G1[j+1,1])
f[nrow(C),1]<-1-sum(f[1:nrow(C)-1,1])
}


Gvi1<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxi<-0
for (k in 1:nrow(C)){
if (C[j,3]< C[k,1]) {auxi<-sum(f[k,])+ auxi}
}
Gvi1[j,]<-auxi
}



F0<-t(J)%*%f
Q0<-Gvi1/F0
for(i in 1:nrow(h)){
h[i,]<- 1/(NJ[i,]+(J[i,])%*%Q0)
}

S0<-max(abs(f-f0))
f0<-f

} 

indbb<-seq(1,nrow(C),by=1)
indbb9<-seq(1,nrow(C),by=1)

kk0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb9<-(C[,1]==C[i,1])
pos9<-min(which(indbb9==TRUE))
if(pos9==1){
kk0b[indbb9]<-sum(f[indbb9])}
if(pos9>1){
kk0b[indbb9]<-sum(f[indbb9])
}
}


mult4 <- tabulate(match(C[,1],unique(C[,1])))
if(sum(mult4)==length(unique(C[,1]))){   
Fval <- (kk0b)}
if(sum(mult4)>length(unique(C[,1]))){
weigth4<-kk0b[!duplicated(C[,1])]
Fval<- (weigth4)}




hh0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb<-(C[,1]==C[i,1])
pos<-min(which(indbb==TRUE))
if(pos==1){
hh0b[indbb]<-sum(h[indbb])}
if(pos>1){
hh0b[indbb]<-sum(h[indbb])
}
}


mult5 <- tabulate(match(C[,1],unique(C[,1])))
if(sum(mult5)==length(unique(C[,1]))){   
hval <- (hh0b)}
if(sum(mult5)>length(unique(C[,1]))){
weigth5<-hh0b[!duplicated(C[,1])]
hval<- (weigth5)}



x<-unique(C[,1])
events<-sum(mult4)
n.event<-mult4
f<-Fval
h<-hval
FF<-cumsum(f)
Sob<-1-FF+f
Sob[Sob<1e-12]<-0



if (boot==TRUE){ 
if (is.na(B)==TRUE) B<-500

if(trunc=="double"|trunc=="left"){


M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
ind<-seq(1,nrow(C),by=1)
indbb<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)


if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for (b in 1:B){

indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
M2b<-matrix(0,nrow=nrow(M1b),ncol=ncol(M1b))
ord<-order(M1b[,1])
M2b[,1]<-sort(M1b[,1])
M2b[,2:ncol(M1b)]<-M1b[ord,2:ncol(M1b)]

NJb<-matrix(data=0,ncol=1,nrow=nrow(C))
Kb<-NJb
for(j in 1:nrow(C)){
for(k in 1:nrow(C)){
if (M2b[k,2]<=M2b[j,1]& M2b[k,1]>=M2b[j,1]) {NJb[j,1]<-NJb[j,1]+1}
}}
Kb[,1]<-rep(1,length=nrow(C))


mNJb<-min(NJb[1:nrow(C)-1,1])
if (mNJb!=1) {Kb[,1]<-NJb[,1]}
if (mNJb==1){
istarb<-max(which(NJb[1:nrow(C)-1,1]==1))
Kb[1:istarb,1]<-2
if(istarb<nrow(M2b)){Kb[(istarb+1):nrow(M2b),1]<-NJb[(istarb+1):nrow(M2b),1]}
}
NJb<-as.matrix(apply(cbind(NJb,Kb),1,max),ncol=1)


Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (k in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[k,1]>=a2 & M2b[k,1]<=b2) Jb[k,j]<-1
}}
h0b<-1/NJb

G0b<-matrix(data=0,ncol=1,nrow=nrow(C))
G0b[1,1]<-1
for(j in 2:nrow(C)){
G0b[j,1]<-exp(sum(log(1-h0b[1:(j-1),1])))
}
f0b<-matrix(data=G0b[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f0[j,1]<-(G0b[j,1]-G0b[j+1,1])
}
Gvi0b<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxib<-0
for (k in 1:nrow(C)){
if (M2b[j,3]< M2b[k,1]) {auxib<-sum(f0b[k,])+ auxib}
}
Gvi0b[j,]<-auxib
}
F0b<-t(Jb)%*%f0b
Q0b<-Gvi0b/F0b
h1b<-h0b
for(i in 1:nrow(h0b)){
h1b[i,]<- 1/(NJb[i,]+(Jb[i,])%*%Q0b)
}

S0b<-1
if(is.na(error)==TRUE) error<-1e-6

if(is.na(nmaxit)==TRUE)nmaxit<-100

iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
G1b<-matrix(data=0,ncol=1,nrow=nrow(C))
G1b[1,1]<-1
for(j in 2:nrow(C)){
G1b[j,1]<-exp(sum(log(1-h1b[1:(j-1),1])))
}
f1b<-matrix(data=G1b[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f1b[j,1]<-(G1b[j,1]-G1b[j+1,1])
f1b[nrow(C),1]<-1-sum(f1b[1:nrow(C)-1,1])
}
Gvi1b<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxib<-0
for (k in 1:nrow(C)){
if (M2b[j,3]< M2b[k,1]) {auxib<-sum(f1b[k,])+ auxib}
}
Gvi1b[j,]<-auxib
}
F0b<-t(Jb)%*%f1b
Q0b<-Gvi1b/F0b
for(i in 1:nrow(h1b)){
h1b[i,]<- 1/(NJb[i,]+(Jb[i,])%*%Q0b)
}
S0b<-max(abs(f1b-f0b))
f0b<-f1b


}

ff0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb1<-(C[,1]==C[i,1])
pos1<-min(which(indbb1==TRUE))
if(pos1==1){
ff0b[indbb1]<-sum(f1b[indbb1])}
if(pos1>1){
ff0b[indbb1]<-sum(f1b[indbb1])
}
}



FF0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb<-(C[,1]==C[i,1])
pos<-min(which(indbb==TRUE))
if(pos==1){
FF0b[indbb]<-sum(f1b[indbb])}
if(pos>1){
FF0b[indbb]<-sum(f1b[indbb])+FF0b[pos-1]
}
}
Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(ff0b)
M_IF0Sob[,b]<-as.vector(Sobb)
}
}


if(trunc=="right"){

M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
ind<-seq(1,nrow(C),by=1)
indbb2<-seq(1,nrow(C),by=1)
indbbb<-seq(1,nrow(C),by=1)

if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for (b in 1:B){


indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
M2b<-matrix(0,nrow=nrow(M1b),ncol=ncol(M1b))
ord<-order(M1b[,1])
M2b[,1]<-sort(M1b[,1])
M2b[,2:ncol(M1b)]<-M1b[ord,2:ncol(M1b)]

NJb<-matrix(data=0,ncol=1,nrow=nrow(C))
Kb<-NJb
for(j in 1:nrow(C)){
for(k in 1:nrow(C)){
if (M2b[k,2]<=M2b[j,1]& M2b[k,1]>=M2b[j,1]) {NJb[j,1]<-NJb[j,1]+1}
}}
Kb[,1]<-rep(1,length=nrow(C))


mNJb<-min(NJb[1:nrow(C)-1,1])
if (mNJb!=1) {Kb[,1]<-NJb[,1]}
if (mNJb==1){
istarb<-max(which(NJb[1:nrow(C)-1,1]==1))
Kb[1:istarb,1]<-2
if(istarb<nrow(M2b)){Kb[(istarb+1):nrow(M2b),1]<-NJb[(istarb+1):nrow(M2b),1]}
}
NJb<-as.matrix(apply(cbind(NJb,Kb),1,max),ncol=1)


Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (k in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[k,1]>=a2 & M2b[k,1]<=b2) Jb[k,j]<-1
}}
h0b<-1/NJb

G0b<-matrix(data=0,ncol=1,nrow=nrow(C))
G0b[1,1]<-1
for(j in 2:nrow(C)){
G0b[j,1]<-exp(sum(log(1-h0b[1:(j-1),1])))
}
f0b<-matrix(data=G0b[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f0[j,1]<-(G0b[j,1]-G0b[j+1,1])
}
Gvi0b<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxib<-0
for (k in 1:nrow(C)){
if (M2b[j,3]< M2b[k,1]) {auxib<-sum(f0b[k,])+ auxib}
}
Gvi0b[j,]<-auxib
}
F0b<-t(Jb)%*%f0b
Q0b<-Gvi0b/F0b
h1b<-h0b
for(i in 1:nrow(h0b)){
h1b[i,]<- 1/(NJb[i,]+(Jb[i,])%*%Q0b)
}

 S0b<-1
if(is.na(error)==TRUE) error<-1e-6

if(is.na(nmaxit)==TRUE)nmaxit<-100

iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
G1b<-matrix(data=0,ncol=1,nrow=nrow(C))
G1b[1,1]<-1
for(j in 2:nrow(C)){
G1b[j,1]<-exp(sum(log(1-h1b[1:(j-1),1])))
}
f1b<-matrix(data=G1b[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f1b[j,1]<-(G1b[j,1]-G1b[j+1,1])
f1b[nrow(C),1]<-1-sum(f1b[1:nrow(C)-1,1])
}
Gvi1b<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxib<-0
for (k in 1:nrow(C)){
if (M2b[j,3]< M2b[k,1]) {auxib<-sum(f1b[k,])+ auxib}
}
Gvi1b[j,]<-auxib
}
F0b<-t(Jb)%*%f1b
Q0b<-Gvi1b/F0b
for(i in 1:nrow(h1b)){
h1b[i,]<- 1/(NJb[i,]+(Jb[i,])%*%Q0b)
}
S0b<-max(abs(f1b-f0b))
f0b<-f1b


}




M2b[,1]<--M2b[,1]
ord2<-order(M2b[,1])
M2b[,1]<-sort(M2b[,1])
E<--C
ord7<-order(E[,1])
E[,1]<-sort(E[,1])

f1b<-f1b[ord2]

FF0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbbb<-(E[,1]==E[i,1])
pos1<-min(which(indbbb==TRUE))
if(pos1==1){
FF0b[indbbb]<-sum(f1b[indbbb])}
if(pos1>1){
FF0b[indbbb]<-sum(f1b[indbbb])+FF0b[pos1-1]
}
}

fb<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb2<-(E[,1]==E[i,1])
pos2<-min(which(indbb2==TRUE))
if(pos2==1){
fb[indbb2]<-sum(f1b[indbb2])}
if(pos2>1){
fb[indbb2]<-sum(f1b[indbb2])
}
}


Sobb<-1-FF0b+fb
Sobb[Sobb<1e-12]<-0
FF0b[FF0b<1e-12]<-0

M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(fb)
M_IF0Sob[,b]<-as.vector(Sobb)

}
}


M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]

} 



if(boot==TRUE){

if(trunc=="double"|trunc=="left"){

if((display.F==TRUE)&(display.S==TRUE)){

dev.new()
par(mfrow=c(1,2))


plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
}


plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
}


}



if((display.S==FALSE)&(display.F==TRUE)){
dev.new()
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(x),1,lty=3)
segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
}
}

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(x),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
}
}}

if(trunc=="right"){

x<--unique(C[,1])
f<-Fval


ord1<-order(x)
x<-sort(x)
mult4<-mult4[ord1]
n.event<-mult4

f<-f[ord1]
FF<-cumsum(f)

FF[FF<1e-12]<-0
Sob<-1-FF+f
h<-f/Sob


if((display.F==TRUE)&(display.S==TRUE)){

dev.new()
par(mfrow=c(1,2))

C[,1]<--C[order(-C[,1])]

plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(min(C[,1]),0,min(C[,1]),upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperF[i],C[,1][i+1],upperF[i], lty=3)
segments(C[,1][i+1],upperF[i],C[,1][i+1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(min(C[,1]),0,min(C[,1]),lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerF[i],C[,1][i+1],lowerF[i], lty=3)
segments(C[,1][i+1],lowerF[i],C[,1][i+1],lowerF[i+1],lty=3)
}


plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))



for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperS[i],C[,1][i+1],upperS[i], lty=3)
segments(C[,1][i+1],upperS[i],C[,1][i+1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerS[i],C[,1][i+1],lowerS[i], lty=3)
segments(C[,1][i+1],lowerS[i],C[,1][i+1],lowerS[i+1],lty=3)

}

}

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(-C[order(-C[,1])])-(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),1,-C[order(-C[,1])][1],1,lty=3)
segments(max(-C[order(-C[,1])]),0,max(-C[order(-C[,1])])+(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),0,lty=3)
segments(max(-C[order(-C[,1])]),0,max(-C[order(-C[,1])]),min(upperS), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(-C[order(-C[,1])][i],upperS[i],-C[order(-C[,1])][i+1],upperS[i], lty=3)
segments(-C[order(-C[,1])][i+1],upperS[i],-C[order(-C[,1])][i+1],upperS[i+1],lty=3)
}


segments(min(-C[order(-C[,1])])-(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),1,-C[order(-C[,1])][1],1,lty=3)
segments(max(-C[order(-C[,1])]),0,max(-C[order(-C[,1])])+(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(x),0,lty=3)
segments(max(-C[order(-C[,1])]),0,max(-C[order(-C[,1])]),min(lowerS), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(-C[order(-C[,1])][i],lowerS[i],-C[order(-C[,1])][i+1],lowerS[i], lty=3)
segments(-C[order(-C[,1])][i+1],lowerS[i],-C[order(-C[,1])][i+1],lowerS[i+1],lty=3)
}
}

if((display.F==TRUE)&(display.S==FALSE)){
dev.new()
par(mfrow=c(1,1))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(-C[order(-C[,1])])-(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),0,-C[order(-C[,1])][1],0,lty=3)
segments(max(-C[order(-C[,1])]),1,max(-C[order(-C[,1])])+(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),1,lty=3)
segments(-C[order(-C[,1])][1],0,-C[order(-C[,1])][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(-C[order(-C[,1])][i],upperF[i],-C[order(-C[,1])][i+1],upperF[i], lty=3)
segments(-C[order(-C[,1])][i+1],upperF[i],-C[order(-C[,1])][i+1],upperF[i+1],lty=3)
}


segments(min(-C[order(-C[,1])])-(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),0,-C[order(-C[,1])][1],0,lty=3)
segments(max(-C[order(-C[,1])]),1,max(-C[order(-C[,1])])+(max(-C[order(-C[,1])])-min(-C[order(-C[,1])]))/length(C[,1]),1,lty=3)
segments(-C[order(-C[,1])][1],0,-C[order(-C[,1])][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(-C[order(-C[,1])][i],lowerF[i],-C[order(-C[,1])][i+1],lowerF[i], lty=3)
segments(-C[order(-C[,1])][i+1],lowerF[i],-C[order(-C[,1])][i+1],lowerF[i+1],lty=3)
}

}

}
}
if (boot==FALSE|B<40){


if(trunc=="double"|trunc=="left"){

if((display.F==TRUE)&(display.S==TRUE)){

dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}




plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")
segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

}



if((display.S==FALSE)&(display.F==TRUE)){
dev.new()
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}

}

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")
segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}



}}

if(trunc=="right"){

x<--unique(C[,1])
events<-sum(mult4)
n.event<-mult4
f<-Fval

ord1<-order(x)
x<-sort(x)
mult4<-mult4[ord1]
n.event<-mult4
f<-f[ord1]
FF<-cumsum(f)
Sob<-1-FF+f
h<-f/Sob
Sob[Sob<1e-12]<-0
FF[FF<1e-12]<-0


if((display.F==TRUE)&(display.S==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}

plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")
segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}



 }

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()

par(mfrow=c(1,1))
plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob))


for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}
 }

if((display.F==TRUE)&(display.S==FALSE)){
dev.new()
par(mfrow=c(1,1))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}
 }

}

}


cat("n.iterations",iter,"\n")
cat("S0",S0,"\n")
cat("events",events,"\n")
if(boot==TRUE){
 cat("B",B,"\n")
 cat("alpha",alpha,"\n")

summary<-cbind("time"=x,"n.event"=mult4,"density"=round(f,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5), "hazard"=round(h,5))

colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival", "hazard")
rownames(summary)<-rep("",times=length(x))
print(summary,digits=5, justify="left")
return(invisible(list(n.iterations=iter, events=events, B=B, alpha=alpha,time=x, n.event=mult4, density=round(as.vector(f),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), hazard=round(as.vector(h),5), NJ=as.vector(NJ), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
lower.Sob=round(lowerS,5))))



}

if(boot==FALSE){
 
summary<-cbind("time"=x,"n.event"=mult4,"density"=round(f,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5), "hazard"=round(h,5))

colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival", "hazard")
rownames(summary)<-rep("",times=length(x))
print(summary,digits=5, justify="left")
return(invisible(list(n.iterations=iter, events=events, time=x, n.event=mult4, density=round(as.vector(f),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5),  hazard=round(as.vector(h),5), NJ=as.vector(NJ))))


}
}

