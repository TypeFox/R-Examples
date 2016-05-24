print.projection<-function(x,...){
    attributes(x)$class<-NULL
    print(x)
}
