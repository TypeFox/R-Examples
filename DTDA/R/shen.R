shen <-
function(X, U=NA, V=NA, wt=NA, error=NA, nmaxit=NA,
 boot=TRUE, boot.type="simple",B=NA, alpha=NA, 
display.FS=FALSE, display.UV=FALSE,
 plot.joint=FALSE, plot.type=NULL){

trunc<-"both"
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

if(trunc=="both"){
if(any(is.na(U))==TRUE|any(is.na(V))==TRUE|any(is.na(X))==TRUE){
navec<-c(which(is.na(X)),which(is.na(U)),which(is.na(V)))
X<-X[-navec]
U<-U[-navec]
V<-V[-navec]
}}



if(is.na(wt)==TRUE) wt<-rep(1/length(X),times=length(X))

D<-cbind(X,U,V)
if(all(is.na(V))==TRUE){
D[,3]<-rep(max(D[,1])+1,length(X))
}

if(all(is.na(U))==TRUE){
D[,2]<-rep(min(D[,1])-1,length(X))
D[,1]<-D[,1]
D[,3]<-D[,3]

}


ord<-order(D[,1])
C<-matrix(0,nrow=nrow(D),ncol=ncol(D))
EE<-matrix(0,nrow=nrow(D),ncol=ncol(D))
C[,1]<-sort(D[,1])
C[,2:ncol(D)]<-D[ord,2:ncol(D)]

if(is.na(error)==TRUE) error<-1e-6
J<-matrix(data=0,ncol=nrow(C),nrow=nrow(C))
for (i in 1:nrow(C)){
for(j in 1:nrow(C)) {
a1<-min(C[j,2],C[j,3])
      b1<-max(C[j,2],C[j,3])
if(C[i,1]>=a1 & C[i,1]<=b1) J[i,j]<-1
}}

JI<-t(J)
f0<-matrix(data=1/nrow(C),ncol=1,nrow=nrow(C))
f<-f0
S0<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iter<-0
while(S0>error | iter>nmaxit){
iter<-iter+1
if (iter>nmaxit) stop("Default number of iterations not enough for convergence")
F0<-JI%*%f
k0<-((sum(1/F0))^(-1))*(1/F0)
if(sum(k0)!=1)k0<-k0/sum(k0)
k<-k0
K0<-J%*%k
f<-((sum(1/K0))^(-1))*(1/K0)
if(sum(f)!=1)f<-f/sum(f)
S0<-max(abs(f-f0))
f0<-f
k0<-k
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

FF<-cumsum(Fval)
Sob<-1-FF+Fval
Sob[Sob<1e-12]<-0

if(trunc=="both"){

indbbb<-seq(1,nrow(C),by=1)
indbbu<-seq(1,nrow(C),by=1)
indbbv<-seq(1,nrow(C),by=1)



kMUV<-cbind(U,V,k)

ordUV<-order(kMUV[,1])
kMUV[,1]<-sort(kMUV[,1])
kMUV[,2]<-kMUV[ordUV,2]


kuv<-numeric(nrow(C))
for(i in 1:nrow(kMUV)){
indbbb<-((kMUV[,1]==kMUV[i,1])&(kMUV[,2]==kMUV[i,2]))
pos1<-min(which(indbbb==TRUE))
if(pos1==1){
kuv[indbbb]<-sum(k[indbbb])}
if(pos1>1){
kuv[indbbb]<-sum(k[indbbb])
}
}

KUUVV<-cbind(kMUV[,1],kMUV[,2])

multC<-dim(unique(KUUVV,margin=1))[1] 
if(multC==dim(KUUVV)[1]){
mult8<-rep(1,times=multC)   
FKval <- (kuv*mult8)}
if(multC<dim(KUUVV)[1]){
weigth8<-kuv[!duplicated(KUUVV,margin=1)]
mult8<-numeric(multC)
for(i in 1:multC){
for(j in 1:dim(KUUVV)[1]){
aux8<-sum(KUUVV[j,]==unique(KUUVV,margin=1)[i,])
if(aux8==2) mult8[i]<-mult8[i]+1
}}
FKval<- (weigth8)}


kMU<-cbind(U,k)
ordU<-order(kMU[,1])
kMU[,1]<-sort(kMU[,1])
kMU[,2]<-kMU[ordU,2]


kk0u<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbbu<-(kMU[,1]==kMU[i,1])
posu<-min(which(indbbu==TRUE))
if(posu==1){
kk0u[indbbu]<-sum(kMU[,2][indbbu])}
if(posu>1){
kk0u[indbbu]<-sum(kMU[,2][indbbu])
}
}


multu <- tabulate(match(kMU[,1],unique(kMU[,1])))
if(sum(multu)==length(unique(kMU[,1]))){   
fUval <- (kk0u)}
if(sum(multu)>length(unique(kMU[,1]))){
weigthu<-kk0u[!duplicated(kMU[,1])]
fUval<- (weigthu)}


UU<-unique(kMU[,1])
fU<-cumsum(fUval)


kMV<-cbind(V,k)
ordV<-order(kMV[,1])
kMV[,1]<-sort(kMV[,1])
kMV[,2]<-kMV[ordV,2]

kk0v<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbbv<-(kMV[,1]==kMV[i,1])
posv<-min(which(indbbv==TRUE))
if(posv==1){
kk0v[indbbv]<-sum(kMV[,2][indbbv])}
if(posv>1){
kk0v[indbbv]<-sum(kMV[,2][indbbv])
}
}


multv <- tabulate(match(kMV[,1],unique(kMV[,1])))
if(sum(multv)==length(unique(kMV[,1]))){   
fVval <- (kk0v)}
if(sum(multv)>length(unique(kMV[,1]))){
weigthv<-kk0v[!duplicated(kMV[,1])]
fVval<- (weigthv)}

VV<-unique(kMV[,1])
fV<-cumsum(fVval)

ordUV<-order(kMUV[,1])
EE[,1]<-X[ordUV]


KKG<-matrix(data=0, ncol=1, nrow=nrow(EE))
for (i in 1:nrow(EE)){
for (j in 1:nrow(EE)){
if(EE[i,1]>=kMUV[j,1]& EE[i,1]<=kMUV[j,2])
KKG[i,]<-KKG[i,]+kMUV[j,3]
}}

multG <- tabulate(match(EE[,1],unique(EE[,1])))
if(sum(multG)==length(unique(EE[,1]))){   
fGval <- (KKG)}
if(sum(multG)>length(unique(EE[,1]))){
weigthG<-KKG[!duplicated(EE[,1])]
fGval<- (weigthG)}


ordc<-order(unique(EE[,1]))
biasf<-fGval[ordc]



KK<-matrix(data=0, ncol=nrow(C), nrow=nrow(C))
for (i in 1:nrow(C)){
for (j in 1:nrow(C)){
for (l in 1:nrow(C)){
if(kMU[i,1]>=kMUV[l,1]& kMV[j,1]>=kMUV[l,2])
KK[i,j]<-KK[i,j]+kMUV[l,3]
}}}


if(plot.joint==TRUE){
if(plot.type=="image"){
sU<-sort(U)
sV<-sort(V)
if(any(diff(sU)==0)& all(diff(sV)!=0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
image(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}

if(all(diff(sU)!=0)& any(diff(sV)==0)){
stepV<-diff(range(V))/(length(V)*length(unique(V)))
image(sort(U),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}

if(any(diff(sU)==0)& any(diff(sV)==0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
stepV<-diff(range(V))/(length(V)*length(unique(V)))
image(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}

if(all(diff(sU)!=0)& all(diff(sV)!=0)){
image(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}}

if(plot.type=="persp"){
fcol<-topo.colors(10)[cut(KK[2:nrow(C),2:nrow(C)],10,include.lowest=TRUE)]
sU<-sort(U)
sV<-sort(V)
if(any(diff(sU)==0)& all(diff(sV)!=0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
if(all(diff(sU)!=0)& any(diff(sV)==0)){
stepV<-diff(range(V))/(length(V)*length(unique(V)))
persp(x=sort(U),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
if(any(diff(sU)==0)& any(diff(sV)==0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
stepV<-diff(range(V))/(length(V)*length(unique(V)))
persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
if(all(diff(sU)!=0)& all(diff(sV)!=0)){
persp(x=sort(U),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
}


}
if(plot.joint==FALSE){
}

}

if(trunc=="left"){

indbbb<-seq(1,nrow(C),by=1)
indbb9<-seq(1,nrow(C),by=1)


kMU<-cbind(U,k)
ordU<-order(kMU[,1])
kMU[,1]<-sort(kMU[,1])
kMU[,2]<-kMU[ordU,2]



EE[,1]<-X[ordU]


KKG<-matrix(data=0, ncol=1, nrow=nrow(EE))
for (i in 1:nrow(EE)){
for (j in 1:nrow(EE)){
if(EE[i,1]>=kMU[j,1])
KKG[i,]<-KKG[i,]+kMU[j,2]
}}

multG <- tabulate(match(EE[,1],unique(EE[,1])))
if(sum(multG)==length(unique(EE[,1]))){   
fGval <- (KKG)}
if(sum(multG)>length(unique(EE[,1]))){
weigthG<-KKG[!duplicated(EE[,1])]
fGval<- (weigthG)}


ordc<-order(unique(EE[,1]))
biasf<-fGval[ordc]


kk0b<-numeric(nrow(kMU))
for(i in 1:nrow(kMU)){
indbb9<-(kMU[,1]==kMU[i,1])
pos9<-min(which(indbb9==TRUE))
if(pos9==1){
kk0b[indbb9]<-sum(kMU[,2][indbb9])}
if(pos9>1){
kk0b[indbb9]<-sum(kMU[,2][indbb9])
}
}



mult3 <- tabulate(match(kMU[,1],unique(kMU[,1])))
if(sum(mult3)==length(unique(kMU[,1]))){   
fUval <- (kk0b)}
if(sum(mult3)>length(unique(kMU[,1]))){
weigth3<-kk0b[!duplicated(kMU[,1])]
fUval<- (weigth3)}

UU<-unique(kMU[,1])
ku<-fUval
fU<-cumsum(fUval)


}


if(trunc=="right"){



indbbb<-seq(1,nrow(C),by=1)
indbb9<-seq(1,nrow(C),by=1)


kMV<-cbind(V,k)
ordV<-order(kMV[,1])
kMV[,1]<-sort(kMV[,1])
kMV[,2]<-kMV[ordV,2]


EE[,1]<-X[ordV]


KKG<-matrix(data=0, ncol=1, nrow=nrow(EE))
for (i in 1:nrow(EE)){
for (j in 1:nrow(EE)){
if(EE[i,1]<=kMV[j,1])
KKG[i,]<-KKG[i,]+kMV[j,2]
}}

multG <- tabulate(match(EE[,1],unique(EE[,1])))
if(sum(multG)==length(unique(EE[,1]))){   
fGval <- (KKG)}
if(sum(multG)>length(unique(EE[,1]))){
weigthG<-KKG[!duplicated(EE[,1])]
fGval<- (weigthG)}


ordc<-order(unique(EE[,1]))
biasf<-fGval[ordc]


kk0b<-numeric(nrow(kMV))
for(i in 1:nrow(kMV)){
indbb9<-(kMV[,1]==kMV[i,1])
pos9<-min(which(indbb9==TRUE))
if(pos9==1){
kk0b[indbb9]<-sum(kMV[,2][indbb9])}
if(pos9>1){
kk0b[indbb9]<-sum(kMV[,2][indbb9])
}
}


mult4 <- tabulate(match(kMV[,1],unique(kMV[,1])))
if(sum(mult4)==length(unique(kMV[,1]))){   
fVval <- (kk0b)}
if(sum(mult4)>length(unique(kMV[,1]))){
weigth4<-kk0b[!duplicated(kMV[,1])]
fVval<- (weigth4)}


VV<-unique(kMV[,1])
kv<-fVval
fV<-cumsum(kv)


}


if(boot==TRUE){

if(boot.type=="simple"){
if (is.na(B)==TRUE) B<-500

if(trunc=="both"){

M1b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)
kMVb<-matrix(0,nrow=nrow(C),ncol=B)
M_fU<-matrix(0,nrow=nrow(C),ncol=B) 
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M3b<-matrix(0,nrow=nrow(C),ncol=B)
M4b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


ind<-seq(1,nrow(C),by=1)
indbb7<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)


if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for (b in 1:B){


indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
ord<-order(M1b[,1])
M2b[,1]<-sort(M1b[,1])
M2b[,2:ncol(M1b)]<-M1b[ord,2:ncol(M1b)]

Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (i in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[i,1]>=a2 & M2b[i,1]<=b2) Jb[i,j]<-1
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
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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
indbb7<-(C[,1]==C[i,1])
pos0<-min(which(indbb7==TRUE))
if(pos0==1){
FF0b[indbb7]<-sum(f1b[indbb7])}
if(pos0>1){
FF0b[indbb7]<-sum(f1b[indbb7])+FF0b[pos0-1]
}
}
Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordA<-order(M2b[,2])
ordB<-order(M2b[,3])
ord3<-sort(M2b[,2])
ord4<-sort(M2b[,3])
M3b[,b]<-ord3
M4b[,b]<-ord4

k1bU[,b]<-k1b[ordU]
k1bV[,b]<-k1b[ordV]

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
 M_fU[j,b]<-sum(k1bU[1:i,b])
 }}

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}

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


M_fU_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]

M_fV_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]

}


if(trunc=="left"){

M1b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)

M_fU<-matrix(0,nrow=nrow(C),ncol=B) 
M3b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 


ind<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)
indbb7<-seq(1,nrow(C),by=1)

if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for (b in 1:B){


indb<-sample(ind,nrow(C),replace=TRUE)


M1b<-C[indb,]
ord<-order(M1b[,1])
M2b[,1]<-sort(M1b[,1])
M2b[,2:ncol(M1b)]<-M1b[ord,2:ncol(M1b)]

Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (i in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[i,1]>=a2 & M2b[i,1]<=b2) Jb[i,j]<-1
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
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
}

ff0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb7<-(C[,1]==C[i,1])
pos7<-min(which(indbb7==TRUE))
if(pos7==1){
ff0b[indbb7]<-sum(f1b[indbb7])}
if(pos7>1){
ff0b[indbb7]<-sum(f1b[indbb7])
}
}



FF0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb1<-(C[,1]==C[i,1])
pos0<-min(which(indbb1==TRUE))
if(pos0==1){
FF0b[indbb1]<-sum(f1b[indbb1])}
if(pos0>1){
FF0b[indbb1]<-sum(f1b[indbb1])+FF0b[pos0-1]
}
}
Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordA<-order(M2b[,2])
ord3<-sort(M2b[,2])
M3b[,b]<-ord3
k1bU[,b]<-k1b[ordU]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}



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


M_fU_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]

}

if(trunc=="right"){

M1b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2bb<-matrix(0,nrow=nrow(C),ncol=ncol(C))

M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMVb<-matrix(0,nrow=nrow(C),ncol=B)
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M4b<-matrix(0,nrow=nrow(C),ncol=B)
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


ind<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)
indbb7<-seq(1,nrow(C),by=1)

if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}


for (b in 1:B){

indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
ord<-order(M1b[,1])
M2b[,1]<-sort(M1b[,1])
M2b[,2:ncol(M1b)]<-M1b[ord,2:ncol(M1b)]

Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (i in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[i,1]>=a2 & M2b[i,1]<=b2) Jb[i,j]<-1
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
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
}

ff0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb7<-(C[,1]==C[i,1])
pos7<-min(which(indbb7==TRUE))
if(pos7==1){
ff0b[indbb7]<-sum(f1b[indbb7])}
if(pos7>1){
ff0b[indbb7]<-sum(f1b[indbb7])
}
}



FF0b<-numeric(nrow(C))
for(i in 1:nrow(C)){
indbb1<-(C[,1]==C[i,1])
pos0<-min(which(indbb1==TRUE))
if(pos0==1){
FF0b[indbb1]<-sum(f1b[indbb1])}
if(pos0>1){
FF0b[indbb1]<-sum(f1b[indbb1])+FF0b[pos0-1]
}
}



Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0
Sobb[Sobb<1e-12]<-0
FF0b[FF0b<1e-12]<-0


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)


ordB<-order(M2b[,3])
ord4<-sort(M2b[,3])
M4b[,b]<-ord4
k1bV[,b]<-k1b[ordV]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
                        M_fV[j,b]<-sum(k1bV[1:i,b])
                        }}


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


M_fV_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]


}

}




if(boot.type=="obvious"){

if (is.na(B)==TRUE) B<-500

if(trunc=="both"){


W1<-as.vector(f)
W2<-as.vector(k)
ind<-1:nrow(C)
if (is.na(B)==TRUE) B<-500
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)
kMVb<-matrix(0,nrow=nrow(C),ncol=B)
M_fU<-matrix(0,nrow=nrow(C),ncol=B) 
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M3b<-matrix(0,nrow=nrow(C),ncol=B)
M4b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 
k1bV<-matrix(0,nrow=nrow(C),ncol=B)
indbb<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)


if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for(b in 1:B){


DB<-matrix(0,nrow=nrow(C),ncol=3)
DB[,1]<-sample(C[,1],size=nrow(C),replace=TRUE,prob=W1)
indUV<-sample(ind,size=nrow(C),replace=TRUE,prob=W2)
DB[,2:3]<-C[indUV,2:3]

log1<-(DB[,1]>=DB[,2])
log2<-(DB[,1]<=DB[,3])
indlog<-(log1*log2==0) 

while(sum(indlog)>0){
DB[indlog,1]<-sample(C[,1],size=sum(indlog),replace=TRUE,prob=W1)
indUV<-sample(ind,size=sum(indlog),replace=TRUE,prob=W2)
DB[indlog,2:3]<-C[indUV,2:3]
log1<-(DB[,1]>=DB[,2])
log2<-(DB[,1]<=DB[,3])
indlog<-(log1*log2==0)
}

ordB<-order(DB[,1])
DBB<-matrix(0,nrow=nrow(C),ncol=3)
DBB[,1]<-sort(DB[,1])
DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

error<-1e-6
Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
for (i in 1:nrow(DBB)){
for(j in 1:nrow(DBB)) {
a2<-min(DBB[j,2],DBB[j,3])
      b2<-max(DBB[j,2],DBB[j,3])
if(DBB[i,1]>=a2 & DBB[i,1]<=b2) Jb[i,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=1/nrow(DBB),ncol=1,nrow=nrow(DBB))
f1b<-f0b
S0b<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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
posb<-min(which(indbb==TRUE))
if(posb==1){
FF0b[indbb]<-sum(f1b[indbb])}
if(posb>1){
FF0b[indbb]<-sum(f1b[indbb])+FF0b[posb-1]
}
}

Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordA<-order(DBB[,2])
ordB<-order(DBB[,3])
ord3<-sort(DBB[,2])
ord4<-sort(DBB[,3])
M3b[,b]<-ord3
M4b[,b]<-ord4
k1bU[,b]<-k1b[ordU]
k1bV[,b]<-k1b[ordV]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}



}

M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


M_fU_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]


M_fV_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]



}


if(trunc=="left"){




ind<-1:nrow(C)
indbb<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)



W1<-as.vector(f)
W2<-as.vector(k)


if (is.na(B)==TRUE) B<-500
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)
M_fU<-matrix(0,nrow=nrow(C),ncol=B) 

M3b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 


if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for(b in 1:B){


DB<-matrix(0,nrow=nrow(C),ncol=3)
DB[,1]<-sample(C[,1],size=nrow(C),replace=TRUE,prob=W1)
indUV<-sample(ind,size=nrow(C),replace=TRUE,prob=W2)
DB[,2:3]<-C[indUV,2:3]

log1<-(DB[,1]>=DB[,2])
log2<-(DB[,1]<=DB[,3])
indlog<-(log1*log2==0) 

while(sum(indlog)>0){
DB[indlog,1]<-sample(C[,1],size=sum(indlog),replace=TRUE,prob=W1)
indUV<-sample(ind,size=sum(indlog),replace=TRUE,prob=W2)
DB[indlog,2:3]<-C[indUV,2:3]
log1<-(DB[,1]>=DB[,2])
log2<-(DB[,1]<=DB[,3])
indlog<-(log1*log2==0)
}

ordB<-order(DB[,1])
DBB<-matrix(0,nrow=nrow(C),ncol=3)
DBB[,1]<-sort(DB[,1])
DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

error<-1e-6
Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
for (i in 1:nrow(DBB)){
for(j in 1:nrow(DBB)) {
a2<-min(DBB[j,2],DBB[j,3])
      b2<-max(DBB[j,2],DBB[j,3])
if(DBB[i,1]>=a2 & DBB[i,1]<=b2) Jb[i,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=1/nrow(DBB),ncol=1,nrow=nrow(DBB))
f1b<-f0b
S0b<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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
posa<-min(which(indbb==TRUE))
if(posa==1){
FF0b[indbb]<-sum(f1b[indbb])}
if(posa>1){
FF0b[indbb]<-sum(f1b[indbb])+FF0b[posa-1]
}
}
Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0

M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordA<-order(DBB[,2])
ord3<-sort(DBB[,2])
M3b[,b]<-ord3
k1bU[,b]<-k1b[ordU]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}




}

M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


M_fU_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]



}



if(trunc=="right"){


ind<-1:nrow(C)

indbb<-seq(1,nrow(C),by=1)
indbb1<-seq(1,nrow(C),by=1)



W1<-as.vector(f)
W2<-as.vector(k)

if (is.na(B)==TRUE) B<-500
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)

kMVb<-matrix(0,nrow=nrow(C),ncol=B)
 
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M4b<-matrix(0,nrow=nrow(C),ncol=B) 
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


if(B<40){
cat("Warning. Number of replicates less than 40","\n")
cat("Confidence bands cannot be computed","\n")
}



for(b in 1:B){


DB<-matrix(0,nrow=nrow(C),ncol=3)
DB[,1]<-sample(C[,1],size=nrow(C),replace=TRUE,prob=W1)
indUV<-sample(ind,size=nrow(C),replace=TRUE,prob=W2)
DB[,2:3]<-C[indUV,2:3]

log1<-(DB[,1]>=DB[,2])
log2<-(DB[,1]<=DB[,3])
indlog<-(log1*log2==0) 

while(sum(indlog)>0){
DB[indlog,1]<-sample(C[,1],size=sum(indlog),replace=TRUE,prob=W1)
indUV<-sample(ind,size=sum(indlog),replace=TRUE,prob=W2)
DB[indlog,2:3]<-C[indUV,2:3]
log1<-(DB[,1]>=DB[,2])
log2<-(DB[,1]<=DB[,3])
indlog<-(log1*log2==0)
}

ordB<-order(DB[,1])
DBB<-matrix(0,nrow=nrow(C),ncol=3)
DBBB<-matrix(0,nrow=nrow(C),ncol=3)

DBB[,1]<-sort(DB[,1])
DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

error<-1e-6
Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
for (k in 1:nrow(DBB)){
for(j in 1:nrow(DBB)) {
a2<-min(DBB[j,2],DBB[j,3])
      b2<-max(DBB[j,2],DBB[j,3])
if(DBB[k,1]>=a2 & DBB[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=1/nrow(DBB),ncol=1,nrow=nrow(DBB))
f1b<-f0b
S0b<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iterb<-0
while(S0b>error | iterb>nmaxit){
iterb<-iterb+1
if (iterb>nmaxit) stop("Default number of iterations not enough for convergence")
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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
posw<-min(which(indbb==TRUE))
if(posw==1){
FF0b[indbb]<-sum(f1b[indbb])}
if(posw>1){
FF0b[indbb]<-sum(f1b[indbb])+FF0b[posw-1]
}
}


Sobb<-1-FF0b+ff0b
Sobb[Sobb<1e-12]<-0
Sobb[Sobb<1e-12]<-0
FF0b[FF0b<1e-12]<-0


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)


ordB<-order(DBB[,3])
ord4<-sort(DBB[,3])

M4b[,b]<-ord4
k1bV[,b]<-k1b[ordV]


for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}



}
M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


M_fV_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]


}

}
}
if (boot==TRUE){

if(trunc=="both"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))

plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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




plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}

CC<-sort(C[,2])

segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
}


segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
}


plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}


CCC<-sort(C[,3])

segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
segments(sort(C[,3])[1],0,sort(C[,3])[1],upperV[1], lty=3)

for(i in 1:(length(C[,3])-1)){
segments(CCC[i],upperV[i],CCC[i+1],upperV[i], lty=3)
segments(CCC[i+1],upperV[i],CCC[i+1],upperV[i+1],lty=3)
}


segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
segments(sort(C[,3])[1],0,sort(C[,3])[1],lowerV[1], lty=3)

for(i in 1:(length(C[,3])-1)){
segments(CCC[i],lowerV[i],CCC[i+1],lowerV[i], lty=3)
segments(CCC[i+1],lowerV[i],CCC[i+1],lowerV[i+1],lty=3)
}

}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}

CC<-sort(C[,2])

segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
}


segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
}


plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}


CCC<-sort(C[,3])

segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
segments(sort(C[,3])[1],0,sort(C[,3])[1],upperV[1], lty=3)

for(i in 1:(length(C[,3])-1)){
segments(CCC[i],upperV[i],CCC[i+1],upperV[i], lty=3)
segments(CCC[i+1],upperV[i],CCC[i+1],upperV[i+1],lty=3)
}


segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
segments(sort(C[,3])[1],0,sort(C[,3])[1],lowerV[1], lty=3)

for(i in 1:(length(C[,3])-1)){
segments(CCC[i],lowerV[i],CCC[i+1],lowerV[i], lty=3)
segments(CCC[i+1],lowerV[i],CCC[i+1],lowerV[i+1],lty=3)
}
}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="left"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
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

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
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




plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}

CC<-sort(C[,2])

segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
}


segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
}


}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}

CC<-sort(C[,2])

segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
}


segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

for(i in 1:(length(C[,2])-1)){
segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
}


}

if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="right"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))

plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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



plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}

segments(min(V)-(max(V)-min(V))/length(V),0,sort(V)[1],0,lty=3)
segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
segments(sort(V)[1],0,sort(V)[1],upperV[1], lty=3)

CCCC<-sort(V)

for(i in 1:(length(V)-1)){
segments(CCCC[i],upperV[i],CCCC[i+1],upperV[i], lty=3)
segments(CCCC[i+1],upperV[i],CCCC[i+1],upperV[i+1],lty=3)
}


segments(min(V)-(max(V)-min(V))/length(V),0,CCCC[1],0,lty=3)
segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
segments(CCCC[1],0,CCCC[1],lowerV[1], lty=3)

for(i in 1:(length(V)-1)){
segments(CCCC[i],lowerV[i],CCCC[i+1],lowerV[i], lty=3)
segments(CCCC[i+1],lowerV[i],CCCC[i+1],lowerV[i+1],lty=3)
}


}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}

segments(min(V)-(max(V)-min(V))/length(V),0,sort(V)[1],0,lty=3)
segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
segments(sort(V)[1],0,sort(V)[1],upperV[1], lty=3)

CCCC<-sort(V)

for(i in 1:(length(V)-1)){
segments(CCCC[i],upperV[i],CCCC[i+1],upperV[i], lty=3)
segments(CCCC[i+1],upperV[i],CCCC[i+1],upperV[i+1],lty=3)
}


segments(min(V)-(max(V)-min(V))/length(V),0,CCCC[1],0,lty=3)
segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
segments(CCCC[1],0,CCCC[1],lowerV[1], lty=3)

for(i in 1:(length(V)-1)){
segments(CCCC[i],lowerV[i],CCCC[i+1],lowerV[i], lty=3)
segments(CCCC[i+1],lowerV[i],CCCC[i+1],lowerV[i+1],lty=3)
}


}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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

if((display.FS==FALSE)&(display.UV==FALSE)){

}
}

}
if (boot==FALSE|B<40){
if(trunc=="both"){

if(plot.joint==TRUE){
if(plot.type=="image"){
sU<-sort(U)
sV<-sort(V)
if(any(diff(sU)==0)& all(diff(sV)!=0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
image(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}

if(all(diff(sU)!=0)& any(diff(sV)==0)){
stepV<-diff(range(V))/(length(V)*length(unique(V)))
image(sort(U),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}

if(any(diff(sU)==0)& any(diff(sV)==0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
stepV<-diff(range(V))/(length(V)*length(unique(V)))
image(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}

if(all(diff(sU)!=0)& all(diff(sV)!=0)){
image(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}}

if(plot.type=="persp"){
fcol<-topo.colors(10)[cut(KK[2:nrow(C),2:nrow(C)],10,include.lowest=TRUE)]
sU<-sort(U)
sV<-sort(V)
if(any(diff(sU)==0)& all(diff(sV)!=0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
if(all(diff(sU)!=0)& any(diff(sV)==0)){
stepV<-diff(range(V))/(length(V)*length(unique(V)))
persp(x=sort(U),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
if(any(diff(sU)==0)& any(diff(sV)==0)){
stepU<-diff(range(U))/(length(U)*length(unique(U)))
stepV<-diff(range(V))/(length(V)*length(unique(V)))
persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
if(all(diff(sU)!=0)& all(diff(sV)!=0)){
persp(x=sort(U),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
}
}


}
if(plot.joint==FALSE){
}




if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}


plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}

 }

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}


plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}

 }

if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="left"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}

 }

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
segments(UU[1],0,UU[1],fU[1])

for(i in 1:(length(UU)-1)){
segments(UU[i],fU[i],UU[i+1],fU[i])
segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
}
 }

if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="right"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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


plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}
 
}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
segments(VV[1],0,VV[1],fV[1])

for(i in 1:(length(VV)-1)){
segments(VV[i],fV[i],VV[i+1],fV[i])
segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
}
}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

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

if((display.FS==FALSE)&(display.UV==FALSE)){

}
}


}

cat("n.iterations",iter,"\n")
cat("S0",S0,"\n")
cat("events",events,"\n")

if(boot==TRUE){

 cat("B",B,"\n")
 cat("alpha",alpha,"\n")
cat("Boot",boot.type,"\n")

summary<-cbind("time"=x,"n.event"=mult,"density"=round(Fval,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5))

colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival")
rownames(summary)<-rep("",times=length(x))
print(summary,digits=5, justify="left")

if(trunc=="both"){
return(invisible(list(n.iterations=iter,events=events, B=B, alpha=alpha,Boot=boot.type,time=x, n.event=mult, density=round(as.vector(Fval),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
lower.Sob=round(lowerS,5), density.joint=round(as.vector(FKval),5),marginal.U=round(fU,5), marginal.V=round(fV,5), upper.fU=round(upperU,5), lower.fU=round(lowerU,5),upper.fV=round(upperV,5), lower.fV=round(lowerV,5),cumulative.joint=round(KK,5) )))

}

if(trunc=="left"){
return(invisible(list(n.iterations=iter,events=events, B=B, alpha=alpha,Boot=boot.type,time=x, n.event=mult, density=round(as.vector(Fval),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
lower.Sob=round(lowerS,5), density.joint=round(as.vector(ku),5),marginal.U=round(fU,5),upper.fU=round(upperU,5), lower.fU=round(lowerU,5))))
}
if(trunc=="right"){

return(invisible(list(n.iterations=iter,events=events, B=B, alpha=alpha,Boot=boot.type,time=x, n.event=mult, density=round(as.vector(Fval),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
lower.Sob=round(lowerS,5), density.joint=round(as.vector(kv),5), marginal.V=round(fV,5), upper.fV=round(upperV,5), lower.fV=round(lowerV,5) )))

}
}

if(boot==FALSE){
 
summary<-cbind("time"=x,"n.event"=mult,"density"=round(as.vector(Fval),5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5))

colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival")
rownames(summary)<-rep("",times=length(x))
print(summary,digits=5, justify="left")


if(trunc=="both"){
return(invisible(list(n.iterations=iter,events=events, time=x, n.event=mult, density=round(as.vector(Fval),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), density.joint=round(as.vector(FKval),5),marginal.U=round(fU,5), marginal.V=round(fV,5), cumulative.joint=round(KK,5) )))

}

if(trunc=="left"){
return(invisible(list(n.iterations=iter,events=events, time=x, n.event=mult, density=round(as.vector(Fval),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5),biasf=round(biasf,5), density.joint=round(as.vector(k),5), marginal.U=round(fU,5))))
}
if(trunc=="right"){

return(invisible(list(n.iterations=iter,events=events, time=x, n.event=mult, density=round(as.vector(Fval),5), cumulative.df=round(FF,5), survival=
round(as.vector(Sob),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), density.joint=round(as.vector(k),5), marginal.V=round(fV,5) )))

}

}
}

