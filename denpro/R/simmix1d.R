simmix1d<-function(n,M,sig,p,seed){
#Simulates a mixture of l normal distributions in R^1,
#
#n is the sample size
#M is l-vector, rows are the means
#sig is l-vector, for l:th mixture variance
#p is l-vector, proportion for each mixture
#
#returns n*d-matrix
#
set.seed(seed) 
l<-length(M)
d<-1
data<-rnorm(n)        #n-vektori valkoista kohinaa 
for (i in 1:n){
   ehto<-runif(1)
   alku<-0
   loppu<-p[1]
   lippu<-0
   for (j in 1:(l-1)){
      if ((alku<=ehto) && (ehto<loppu)){
         data[i]<-sig[j]*data[i]+M[j]
         lippu<-1
      }
      alku<-alku+p[j]
      loppu<-loppu+p[j+1]
   }      
   if (lippu==0) data[i]<-sig[l]*data[i]+M[l]
}
data<-t(t(data))  #we make a n*1 matrix
return(data)
}
