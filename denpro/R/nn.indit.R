nn.indit<-function(dendat)
{
n<-dim(dendat)[1]
maxk<-n-1
indmat<-matrix(0,n,maxk)

eta<-dist(dendat)
#i<j eta[n*(i-1) - i*(i-1)/2 + j-i]

for (i in 2:(n-1)){
   i1<-seq(1,i-1)
   j1<-i
   irow1<-eta[n*(i1-1) - i1*(i1-1)/2 + j1-i1]
   j2<-seq(i+1,n)
   irow2<-eta[n*(i-1) - i*(i-1)/2 + j2-i]
   irow<-c(irow1,irow2)
   or<-order(irow)
   poisi<-c(seq(1,i-1),seq(i+1,n))
   indmat[i,]<-poisi[or]
}

i<-1
j<-seq(i+1,n)
irow<-eta[n*(i-1) - i*(i-1)/2 + j-i]
or<-order(irow)
poisi<-seq(2,n)
indmat[i,]<-poisi[or]

i<-n
i1<-seq(1,n-1)
j<-i
irow<-eta[n*(i1-1) - i1*(i1-1)/2 + j-i1]
or<-order(irow)
poisi<-seq(1,n-1)
indmat[i,]<-poisi[or]

return(indmat)
}
