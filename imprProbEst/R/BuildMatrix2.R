`BuildMatrix2` <-
function(n,m){
   # CAUTION: This 'n' is rather 'nr'
   A1 <- cbind(diag(n),matrix(0,nrow=n,ncol=m),-diag(n))
   A2 <- cbind(matrix(rep(1,n+m),nrow=1),matrix(rep(0,n),nrow=1))
   A3 <- -A2
   A4 <- -diag(n+m+n)

   rbind(A1,A2,A3,A4)
 }

