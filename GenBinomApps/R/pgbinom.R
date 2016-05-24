pgbinom <-
function(q,size,prob,lower.tail = TRUE,log.p=FALSE){
q<-as.matrix(q)
z=apply(q,1,function(q) sum(dgbinom(c(0:q),size,prob)))
if (lower.tail == TRUE) 
{pgb<-z}
else if (lower.tail == FALSE)
{pgb<-1-z}
if (log.p==TRUE){
pgb<-log(pgb)}
pgb
}
