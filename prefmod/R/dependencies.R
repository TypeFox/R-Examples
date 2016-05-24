# dependencies
dependencies<-function(nobj,Y)
{
     m<-matrix(,nrow=nobj-1,ncol=nobj)

     label.comps<-NULL
     for (j in 2:nobj){
         for (i in 1:(j-1)){
             label.comps<-c(label.comps,paste(i,j,sep=""))
         }
     }

     k<-1
     for (j in 2:nobj){
       for (i in 1:(j-1)){
          m[j-1,i]<-k    #lower triangle
          m[i,j]<-k      #upper triangle
          k<-k+1
       }
     }

     label.intpars<-NULL
     dncol<-(nobj)*(nobj-1)*(nobj-2)/2
     d<-matrix(,ncol=dncol,nrow=nrow(Y))
     l<-1
     for (k in 1:(nobj)){
       for (i in 1:(nobj-2)){
          for (j in (i+1):(nobj-1)){
             d[,l]<-Y[,m[i,k]]*Y[,m[j,k]]
             label.intpars<-c(label.intpars,paste("I",label.comps[m[i,k]],".",label.comps[m[j,k]],sep=""))
             l<-l+1
          }
       }
     }
     list(d=d,label.intpars=label.intpars)
}
