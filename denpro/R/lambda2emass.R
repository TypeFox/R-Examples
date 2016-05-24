lambda2emass<-function(lambda,m,M,sig,p,support=NULL,seed=1,mul=2)
{
#m is the number of Monte Carlo samples
#M is l*d-matrix, rows are the means
#sig is l*d-matrix, for l:th mixture d covariances
#p is l-vector, proportion for each mixture

set.seed(seed)
l<-dim(M)[1]
d<-dim(M)[2]
if (is.null(support)){
   support<-matrix(0,2*d,1)
   for (i in 1:d){
       support[2*i-1]<-min(M[,i]-mul*sig[,i])
       support[2*i]<-max(M[,i]+mul*sig[,i])
   }
}

maksi<-0
for (i in 1:l){
    zig<-sig[i,]
    maksi<-maksi+p[i]*evanor(0)/prod(zig)
}

boxvol<-1
for (i in 1:d) boxvol<-boxvol*(support[2*i]-support[2*i-1])
boxvol<-boxvol*maksi

inside<-0
for (i in 1:m){
    x<-matrix(0,d,1)
    ran<-runif(d+1)
    for (j in 1:d){
        beg<-support[2*j-1]
        end<-support[2*j] 
        x[j]<-beg+(end-beg)*ran[j]
    }
    y<-0+(maksi-0)*ran[d+1]

    arvo<-0
    for (j in 1:l){
        zig<-sig[j,]
        mu<-M[j,]
        arvo<-arvo+p[j]*evanor((x-mu)/zig)/prod(zig)
    }
    
    if ((y<=arvo)&&(y>=lambda)) inside<-inside+1
}

emass<-boxvol*inside/m
return(emass)
}


