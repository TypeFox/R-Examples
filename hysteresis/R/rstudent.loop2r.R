rstudent.loop2r <- function(model,...) {
  g <- model
  n <- length(g$x)
  wr1<- g$x-g$pred.x
  wr2<- g$y-g$pred.y
  var.x <- crossprod(wr1-mean(wr1))/(n-3)
  var.y <- crossprod(wr2-mean(wr2))/(n-3)
  Xmat <- cbind(rep(1,n),sin(g$period.time),cos(g$period.time))
  hX <- Xmat%*%solve(crossprod(Xmat))%*%t(Xmat)
  Ind <- (g$period.time > 0) & (g$period.time < pi)
  if (model$classical==FALSE)
  Ymat <- cbind(rep(1,n),sin(g$period.time)^g$values["m"],cos(g$period.time)^g$values["n"],Ind*sin(g$period.time)^g$values["m"])
  else
    direc <- sign(cos(g$period.time))  
  Ymat <- cbind(rep(1,n),sin(g$period.time)^g$values["m"],direc*abs(cos(g$period.time))^g$values["n"],Ind*sin(g$period.time)^g$values["m"])
  hY <- Ymat%*%solve(crossprod(Ymat))%*%t(Ymat)
  r.Ta <- wr1/sqrt(1-diag(hX))
  r.Tb <- wr2/sqrt(1-diag(hY))     
  wr1 <- r.Ta-mean(r.Ta)
  wr2 <- r.Tb-mean(r.Tb) 
  data.frame("input"=wr1,"output"=wr2)
}
