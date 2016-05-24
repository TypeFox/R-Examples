krylov<-function(A,b,m){
    K<-matrix(,length(b),m)
    dummy<-b
    for (i in 1:m){
        K[,i]<-dummy
        dummy<-A%*%dummy
    }
    return(K)

}