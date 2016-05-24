rstudent.loopsummary <- function(model,...) {
  g <- model
  n <- length(g$x)
  wr1<- g$x-g$pred.x
  wr2<- g$y-g$pred.y
  var.x <- crossprod(wr1-mean(wr1))/(n-3)
  var.y <- crossprod(wr2-mean(wr2))/(n-3)
  Xmat <- cbind(rep(1,n),sin(g$period.time),cos(g$period.time))
  hX <- Xmat%*%solve(crossprod(Xmat))%*%t(Xmat)
  Ymat <- cbind(rep(1,n),sin(g$period.time)^g$Boot.Estimates["m"],cos(g$period.time)^g$Boot.Estimates["n"])
  hY <- Ymat%*%solve(crossprod(Ymat))%*%t(Ymat)
  r.Ta <- wr1/sqrt((1-diag(hX))*var.x)
  r.Tb <- wr2/sqrt((1-diag(hY))*var.y) 
  wr1 <- r.Ta-mean(r.Ta)
  wr2 <- r.Tb-mean(r.Tb) 
  data.frame("input"=wr1,"output"=wr2)
}
