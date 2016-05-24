cal.int <- function(len.b,q,help.env,knots.val) {
   for(i in 1:len.b) {
     INT <- 0
     q.sec <- 1:(q+1)
     for(j in 1:q) {
       y2 <- knots.val$help[i+j]
       y1 <- knots.val$help[i+j-1]
       coef <- get(paste("coef",i,".",j,sep=""),envir=help.env)
       y2 <- 1/q.sec*y2^q.sec
       y1 <- 1/q.sec*y1^q.sec
       INT <- sum(coef*y2)-sum(coef*y1)
       assign(paste("INT",i,".",j,sep=""),INT,envir=help.env)
     }
   }
  return(help.env)
}
