efron.petrosian <-
function(X, U=NA, V=NA, wt=NA, error=NA,
 nmaxit=NA , boot=TRUE, B=NA, alpha=NA, 
display.F=FALSE, display.S=FALSE){

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

if(is.na(wt)==TRUE) wt<-rep(1/length(X),times=length(X))
D<-cbind(X,U,V)

if(all(is.na(U))==TRUE){
D[,1]<-X
D[,2]<-min(D[,1])-1
D[,3]<-V
}
if(all(is.na(V))==TRUE){
D[,3]<-rep(max(D[,1])+1,length(X))
}

ord<-order(D[,1])
C<-matrix(0,nrow=nrow(D),ncol=ncol(D))
C[,1]<-sort(D[,1])
C[,2:ncol(D)]<-D[ord,2:ncol(D)]


if(is.na(error)==TRUE) error<-1e-6
J<-matrix(data=0,ncol=nrow(C),nrow=nrow(C))
for (k in 1:nrow(C)){
for(j in 1:nrow(C)) {
a1<-min(C[j,2],C[j,3])
     b1<-max(C[j,2],C[j,3])
if(C[k,1]>=a1 & C[k,1]<=b1) J[k,j]<-1
}}

JI<-t(J)
f0<-matrix(data=wt,ncol=1,nrow=nrow(C))
f<-f0
S0<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iter<-0
while(S0>error | iter>nmaxit){
iter<-iter+1
if (iter>nmaxit) stop("Default number of iterations not enough for convergence")

F0<-JI%*%f
IF0<-1/F0
If1<-J%*%IF0
f<-1/If1
if(sum(f)!=1)f<-f/sum(f)

S0<-max(abs(f-f0))
f0<-f
}

F0<-JI%*%f
mult <- tabulate(match(C[,1],unique(C[,1])))
if(sum(mult)==length(unique(C[,1]))){   
Fval <- (f*mult)}
if(sum(mult)>length(unique(C[,1]))){
weigth<-f[!duplicated(C[,1])]
Fval<- (weigth*mult)}

x<-unique(C[,1])
events<-sum(mult)
n.event<-mult
f<-Fval
FF<-cumsum(f)
Sob<-1-FF+f
Sob[Sob<1e-12]<-0


if (boot==TRUE){
if (is.na(B)==TRUE) B<-500

if(trunc=="double"|trunc=="left"){


M_IF0<-matrix(0,ncol=B,nrow=nrow(C))
M_IF01<-matrix(0,ncol=B,nrow=nrow(C))
M_IF0Sob<-matrix(0,ncol=B,nrow=nrow(C))
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

Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (k in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[k,1]>=a2 & M2b[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=wt,ncol=1,nrow=nrow(M2b))
f1b<-f0b
 S0b<-1
iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
F0b<-JIb%*%f1b
IF0b<-1/F0b
If1b<-Jb%*%IF0b
f1b<-1/If1b
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
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

M_IF0<-matrix(0,ncol=B,nrow=nrow(C))
M_IF01<-matrix(0,ncol=B,nrow=nrow(C))
M_IF0Sob<-matrix(0,ncol=B,nrow=nrow(C))
ind<-seq(1,nrow(C),by=1)
indbb<-seq(1,nrow(C),by=1)
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

Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (k in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[k,1]>=a2 & M2b[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=wt,ncol=1,nrow=nrow(M2b))
f1b<-f0b
 S0b<-1
iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
F0b<-JIb%*%f1b
IF0b<-1/F0b
If1b<-Jb%*%IF0b
f1b<-1/If1b
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
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


fb<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbbb<-(C[,1]==C[i,1])
pos1<-min(which(indbbb==TRUE))
if(pos1==1){
fb[indbbb]<-sum(f1b[indbbb])}
if(pos1>1){
fb[indbbb]<-sum(f1b[indbbb])
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


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
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

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
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


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
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

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
}
}}

if(trunc=="right"){


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


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(max(C[,1]),1,max(C[,1]),max(upperF), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperF[i],C[,1][i+1],upperF[i], lty=3)
segments(C[,1][i+1],upperF[i],C[,1][i+1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(max(C[,1]),1,max(C[,1]),max(lowerF), lty=3)

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


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperF[i],C[,1][i+1],upperF[i], lty=3)
segments(C[,1][i+1],upperF[i],C[,1][i+1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerF[i],C[,1][i+1],lowerF[i], lty=3)
segments(C[,1][i+1],lowerF[i],C[,1][i+1],lowerF[i+1],lty=3)
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

summary<-cbind("time"=x,"n.event"=mult,"density"=round(f,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5))

colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival")
rownames(summary)<-rep("",times=length(x))
print(summary,digits=5, justify="left")
return(invisible(list(n.iterations=iter, events=events, B=B, alpha=alpha,time=x, n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
lower.Sob=round(lowerS,5))))



}

if(boot==FALSE){
 
summary<-cbind("time"=x,"n.event"=mult,"density"=round(f,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5))

colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival")
rownames(summary)<-rep("",times=length(x))
print(summary,digits=5, justify="left")
return(invisible(list(n.iterations=iter, events=events, time=x, n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5))))


}
}

