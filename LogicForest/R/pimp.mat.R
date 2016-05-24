pimp.mat <-
function(pimps.out, testdata)
{
 tmp.mat<-pimps.out$tmp.mat 
 zero.ids<-c()
 for(i in 1:ncol(tmp.mat))
   {
   ids<-if(all(tmp.mat[,i]==0)) {ids<-i}
   zero.ids<-append(zero.ids, ids)
   }
 if (length(zero.ids) > 0) {tmp.mat<-tmp.mat[,-zero.ids]}
 pimp.ids<-pimps.out$vec.pimpvars 
 subdata<-as.matrix(testdata[,pimp.ids])
 if (is.null(dim(tmp.mat))) {tmp.mat<-matrix(1,1,1)}
 if (nrow(tmp.mat)!=length(pimps.out$vec.primes)) 
   {tmp.mat<-t(tmp.mat)}
 if (is.matrix(tmp.mat)) {npimps<-nrow(tmp.mat)}
 if (is.vector(tmp.mat)) {npimps<-1}
 n<-nrow(subdata)
 pimp.datamat<-matrix(0, nrow=n, ncol=npimps)
 colnames(pimp.datamat)<-pimps.out$vec.primes
 for (i in 1:npimps)
   {
   if (is.matrix(tmp.mat)) {match.matrix<-matrix(0, nrow=n, ncol=ncol(tmp.mat))}
   if (is.vector(tmp.mat)) {match.matrix<-matrix(0, nrow=n, ncol=length(tmp.mat))}
   for (j in 1:n)
     {
     if (is.matrix(tmp.mat)) 
       {
       for (k in 1:ncol(tmp.mat))
         {
         if (tmp.mat[i,k]==1 & subdata[j,k]==1) {match.matrix[j,k]<-1}
         if (tmp.mat[i,k]==-1 & subdata[j,k]==0) {match.matrix[j,k]<-1}
         if (tmp.mat[i,k]==0 & subdata[j,k]==1|tmp.mat[i,k]==0 & subdata[j,k]==0) {match.matrix[j,k]<-1}
         }
       }
    if(is.vector(tmp.mat)) 
       {
       for (k in 1:length(tmp.mat))
         {
         if (tmp.mat[k]==1 & subdata[j,k]==1) {match.matrix[j,k]<-1}
         if (tmp.mat[k]==-1 & subdata[j,k]==0) {match.matrix[j,k]<-1}
         if (tmp.mat[k]==0 & subdata[j,k]==1|tmp.mat[k]==0 & subdata[j,k]==0) {match.matrix[j,k]<-1}
         }
       }
    pimp.datamat[j,i]<-ifelse(all(match.matrix[j,]==1), 1, 0)
    }
   }
 pimp.names<-pimps.out$vec.primes
 pimp.info<-list(pimp.names=pimp.names, pimp.datamat=pimp.datamat)
 pimp.info 
}
