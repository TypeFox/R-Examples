overallmean <-
function(x){
   overall <- mean(x)
   n1 <- dim(x)[1]
   n2 <- dim(x)[2]
   matrix(nrow=n1,ncol=n2,data=rep(overall,n1*n2))
}
