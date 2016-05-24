rstudent.ellipsesummary <- function(model,...) {
  g <- model
  n <- length(g$x)
  wr1<- g$x-g$pred.x
  wr2<- g$y-g$pred.y
  Xmat <- cbind(rep(1,n),sin(g$period.time),cos(g$period.time))
  var.x <- crossprod(wr1-mean(wr1))/(n-3)
  var.y <- crossprod(wr2-mean(wr2))/(n-3)
  h <- Xmat%*%solve(crossprod(Xmat))%*%t(Xmat)
  r.Ta <- wr1/sqrt((1-diag(h))*var.x)
  r.Tb <- wr2/sqrt((1-diag(h))*var.y)  
  wr1 <- r.Ta-mean(r.Ta)
  wr2 <- r.Tb-mean(r.Tb) 
  data.frame("input"=wr1,"output"=wr2)
}
