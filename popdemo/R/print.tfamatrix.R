print.tfamatrix<-function(x,...){
    attributes(x)$class<-NULL
    l<-length(x)
    print(x[1:(l-1)])
}
