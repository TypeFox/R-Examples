montecarlo.ball<-function(dendat,rho,M,seed=1,type="ball")
{
# dendat on n*d matriisi
n<-dim(dendat)[1]
d<-dim(dendat)[2]

if (type=="ball"){
   keski<-colMeans(dendat)
   etais<-matrix(0,n,1)
   for (i in 1:n) etais[i]<-sqrt(sum((keski-dendat[i,])^2))
   masi<-max(etais)
   sade<-masi+rho
   set.seed(seed)
   polap<-sade*sqrt(runif(M))
   polax<-matrix(rnorm(M*d),M,d)
   varia<-matrix(0,M,d)
   for (i in 1:M) varia[i,]<-keski+polap[i]*polax[i,]/sqrt(sum(polax[i,]^2))
}
else {   # type=rectangular
  lows<-matrix(0,d,1)
  higs<-matrix(0,d,1)
  for (i in 1:d){
      ma<-max(dendat[,i])
      mi<-min(dendat[,i])
      lows[i]<-mi-rho
      higs[i]<-ma+rho
  }
  set.seed(seed)
  varia<-matrix(runif(M*d),M,d)
  for (i in 1:d) varia[,i]<-varia[,i]*(higs[i]-lows[i])+lows[i]
}

count<-matrix(0,M,1)
for (i in 1:M){
    point<-varia[i,]
    sisalla<-0
    j<-1
    while ((j<=n)&&(sisalla==0)){
        dista2<-sum((point-dendat[j,])^2)
        if (dista2<=rho^2){ 
            sisalla<-1
            count[i]<-1
        }
        j<-j+1
    }
}

if (type=="ball"){
  voluball<-pi*sade^2
  volu<-voluball*sum(count)/M
}
else{
  volurec<-1
  for (i in 1:d) volurec<-volurec*(higs[i]-lows[i])
  volu<-volurec*sum(count)/M
}

return(volu)
}


