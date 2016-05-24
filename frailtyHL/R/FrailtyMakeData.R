FrailtyMakeData <-
function(y,x,del,z) {
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    zz<-z[[1]]
    if (nrand>1) {
        index1<-nrand
        for (i in 2:index1) zz<-cbind(zz, z[[i]])
    }
    qcum <- cumsum(c(0, q))
    zzz<-matrix(0,n,qcum[nrand+1])
    zzz[1:n,1:qcum[nrand+1]]<-zz[1:n,1:qcum[nrand+1]]
    sort_data<-cbind(y,x,del,zzz)
    sort.res<-sort_data[order(sort_data[,1],na.last=NA),]
    y[1:n,1]<-sort.res[1:n,1]
    index1<-p+1
    x[1:n,1:p]<-sort.res[1:n,2:index1]
    index1<-index1+1
    del[1:n,1]<-sort.res[1:n,index1]
    for (i in 1:nrand) {
        index3<-p+2+qcum[i]+1
        index4<-p+2+qcum[i+1]
        z[[i]][1:n,1:q[i]]<-sort.res[1:n,index3:index4]
    }

t<-matrix(0,n,1)
xx<-matrix(0,n,p)
di<-matrix(0,n,1)
idx1<-0

for (i in 1:n) {
if (del[i,1]==1) {
   idx1<-idx1+1
   t[idx1,1]<-y[i,1]
   for (j in 1:p) {
    xx[idx1,j]<-x[i,j]
   }
}
}

t1<-t
for (i in 1:idx1) {
for (j in 1:idx1) {
if (t1[i,1]==t1[j,1]) {
   if(i != j) {
    t1[j,1]<-0
   }
}
}
}

t2<-matrix(0,idx1,1)
idx2<-0
for (i in 1:idx1) {
if (t1[i,1] != 0) {
 idx2<-idx2+1
 t2[idx2,1]<-t1[i,1]
}
}

di<-matrix(0,idx2,1)
si<-matrix(0,idx2,p)
for (i in 1:idx2) {
   di[i,1]<-0
   for (j in 1:idx1) {
    if(t2[i,1]==t[j,1]) {
     di[i,1]<-di[i,1]+1
     for(k in 1:p) {
      si[i,k]<-si[i,k]+xx[j,k]
     }
    }
   }
}

# Triangle Mat form
Mi<-matrix(0,n,idx2)
for (i in 1:n) {
  t0<-y[i,1]
 for (j in 1:idx2) {
 if (t2[j,1] <= t0) {
  Mi[i,j]=1
 } else Mi[i,j]=0
 }
}
  res<-list(y,x,del,z,Mi,idx2,t2, di)
  return(res)
}

