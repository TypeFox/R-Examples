plot.totpop <-
function(x,...){

hist(x$TOT,xlab="Total Population",probability=TRUE,ylab="Posterior probability",main="",...)
}
