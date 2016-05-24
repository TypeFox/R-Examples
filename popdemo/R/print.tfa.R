print.tfa<-function(x,...){
    attributes(x)$class<-NULL
    print(x)
}
