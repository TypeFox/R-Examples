library(pcalg)

set.seed(43)
nreps <- 50
res <- matrix(FALSE, 3, nreps)
resF <- rep(FALSE, nreps)

for (i in 1:nreps) {
  if (i%%100 == 0) cat("i=",i,"\n")
  set.seed(i)
  p <- 10
  ## generate and draw random DAG :
  myDAG <- randomDAG(p, prob = 0.4)
  myCPDAG <- dag2cpdag(myDAG)
  mcov <- trueCov(myDAG)

  x <- sample(1:10,1)
  y1 <- sample(setdiff(1:10,x),1)
  y2 <- sample(setdiff(1:10,c(x,y1)),1)
  y3 <- sample(setdiff(1:10,c(x,y1,y2)),1)
  ## plot(myCPDAG)
  eff.true1 <- round(causalEffect(myDAG, y1, x),14)
  eff.true2 <- round(causalEffect(myDAG, y2, x),14)
  eff.true3 <- round(causalEffect(myDAG, y3, x),14)
  ## cat("x=",x," y1=",y1," eff=",eff.true1,"\n")
  ## cat("x=",x," y1=",y2," eff=",eff.true2,"\n")

  (eff.est1 <- round(ida(x,y1,mcov,myCPDAG,method="local",verbose=FALSE),14))
  (eff.est2 <- round(ida(x,y2,mcov,myCPDAG,method="local",verbose=FALSE),14))
  (eff.est3 <- round(ida(x,y3,mcov,myCPDAG,method="local",verbose=FALSE),14))
  (eff.estF <- round(idaFast(x,c(y1,y2,y3),mcov,myCPDAG),14))

  res[1,i] <- (eff.true1 %in% eff.est1)
  res[2,i] <- (eff.true2 %in% eff.est2)
  res[3,i] <- (eff.true3 %in% eff.est3)

  resF[i] <- all(eff.estF == rbind(eff.est1,eff.est2,eff.est3))
}

if (!all(res)) stop("Test idaFast: True effects were not recovered!")
if (!all(resF)) stop("Test idaFast: ida and idaFast didn't give same results!")
