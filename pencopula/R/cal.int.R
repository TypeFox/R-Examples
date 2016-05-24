 cal.int <- function(len.k,len.b,q,help.env,knots.val) {
  INT <- matrix(0,q,len.b)
  INT.help <- rep(0,len.b)

  for(i in 1:(len.k-(q-1))) {
    count <- 0
    for(j in 1:q) {
      y2 <- knots.val$val[i+1]
      y1 <- knots.val$val[i]
      coef <- get(paste("coef",i,".",j,sep=""),envir=help.env)
      y2 <- 1/(1:(q+1))*y2^(1:(q+1))
      y1 <- 1/(1:(q+1))*y1^(1:(q+1))
      INT[j,i+count] <- sum(coef*y2)-sum(coef*y1)
      assign(paste("INT",i,".",j,sep=""),INT[j,i+count],envir=help.env)
      count <- count+1
    }
  }
  assign("stand.num",1/colSums(INT),help.env)
}
