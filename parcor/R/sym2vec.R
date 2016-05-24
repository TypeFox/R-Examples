# converts to upper triangle of a symmetric matrix to a vector, probably, this is a very inefficient implementation
sym2vec<-function(A){
    p<-ncol(A)
    v<-c()
    for (i in 1:(p-1)){
        v<-c(v,A[i,(i+1):p])
    }
return(v)


}
