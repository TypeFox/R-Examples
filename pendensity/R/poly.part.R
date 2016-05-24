poly.part <- function(i,j,knots.val,help.env,q,yi=NULL,poly=FALSE) {
  q.sec <- 1:(q+1)
  if(!poly) {
    if(is.null(yi)) return(get(paste("INT",i,".",j,sep=""),envir=help.env))
    else {
      y2 <- yi
      y1 <- knots.val$help[i+j-1]
      coef <- get(paste("coef",i,".",j,sep=""),envir=help.env)
      y2 <- 1/(q.sec)*y2^q.sec
      y1 <- 1/(q.sec)*y1^q.sec
      INT <- sum(coef*y2)-sum(coef*y1)
      return(INT)
    }
  }
  else {
    coef <- get(paste("coef",i,".",j,sep=""),envir=help.env)
    y1 <- knots.val$help[i+j-1]

    y1 <- 1/q.sec*y1^q.sec
    help <- -sum(coef*y1)

    term<-paste("x",q.sec,sep="^")
    term<-paste(1/q.sec,term,sep=" * ")
    term<-paste(coef,term,sep=" * ")
    term <- c(term,help)
    term<-paste(term,collapse=" + ")
    return(term)
  }
} 
