pprodnormalMC<-
  function(q, mu.x, mu.y, se.x, se.y, rho=0, lower.tail=TRUE, n.mc=1e7){
    mu.x <- mu.x/se.x
    mu.y <- mu.y/se.y
    q <- q/(se.x*se.y)
    mean.v <- c(mu.x,mu.y)
    var.mat <- matrix(c(1,rho,rho,1),2)
    a_b <- matrix(rnorm(2*n.mc),ncol=n.mc)
    a_b <- crossprod(chol(var.mat),a_b)+mean.v
    a_b <- t(a_b)
    ab <- a_b[,1]*a_b[,2]
    x <- ab<q
    percentile <- sum(x)/n.mc
    error <- sd(x)/n.mc
    if (!lower.tail)
      percentile <- 1-percentile
    return(list(p=percentile,error=error))
  }
