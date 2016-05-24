plotmtx <-function(m){

tm<-t(m)

nc<-ncol(tm)

m.out<-tm[,(nc:1)]

return(m.out)

}

