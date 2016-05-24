my.bspline <- function(y,K,q,margin.normal=FALSE) {
  knots <- seq(0,1,length=K) 
  if(margin.normal) {
    knots <- qnorm(knots)
    knots[1]<-qnorm(0.0000001)
    knots[length(knots)]<-qnorm(1-0.0000001)
  }
  len.k <- length(knots)
  base.den <- bsplineS(y,breaks=knots,norder=q)
  len.b <- dim(base.den)[2]

  knots.val <- list()
  knots.val$val <- knots
  
  #integration
  help.env <- new.env()
  assign("base.den",base.den,help.env)
  assign("knots.val",knots.val,help.env)
  assign("y",y,help.env)
  assign("q",q,help.env)
  int.my.bspline(help.env)
  stand.num <- get("stand.num",help.env)
  INT <- get("INT",help.env)

  for(j in 1:len.b) base.den[,j] <- base.den[,j]*stand.num[j]

  return(list(base.den=base.den,stand.num=stand.num,knots.val=knots.val,K=K,INT=INT))
}
