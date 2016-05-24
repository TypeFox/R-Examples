vvtz <-
function(v,z){
    dummy<-v%*%(t(v)%*%z)
    dummy<-as.vector(dummy)
    return(dummy)
}
