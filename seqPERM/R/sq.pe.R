sq.pe <-
function(r1,r2,v){
S<-seq(r1,r2,1)
L<-length(S)
b<-lapply(1:v, function(i) c(make(r1,r2,i,v)))
bb<-do.call(cbind,b)
return(bb)}
