createb <-
function(bA,A,bm,m,bM,M){
if(!is.null(A)){
lb<-length(bA)
b<-array(0,c(lb,1))
for(i in 1:lb)
b[i]<-bA[i]
if(!is.null(m)){
lbm<-length(bm)
la<-lb+lbm
ab<-array(0,c(la,1))
for(i in 1:lb)
ab[i]<-bA[i]
for(i in (lb+1):la)
ab[i]<-bm[i-lb]
b<-ab
if(!is.null(M)){
lbM<-length(bM)
lc<-la+lbM
aux<-array(0,c(lc,1))
for(i in 1:la)
aux[i]<-b[i]
for(i in (la+1):lc)
aux[i]<-bM[i-la]
b<-aux
}
}
else if(!is.null(M)){
lbM<-length(bM)
lc<-lb+lbM
aux<-array(0,c(lc,1))
for(i in 1:lb)
aux[i]<-b[i]
for(i in (lb+1):lc)
aux[i]<-bM[i-lb]
b<-aux
}
}
else if(!is.null(m)){
lb<-length(bm)
b<-array(0,c(lb,1))
for(i in 1:lb)
b[i]<-bm[i]
if(!is.null(M)){
lbM<-length(bM)
la<-lb+lbM
aux<-array(0,c(la,1))
for(i in 1:lb)
aux[i]<-b[i]
for(i in (lb+1):la)
aux[i]<-bM[i-lb]
b<-aux
}
}
else if(!is.null(M)){
lb<-length(bM)
b<-array(0,c(lb,1))
for(i in 1:lb)
b[i]<-bM[i]
}
return(b)
}
