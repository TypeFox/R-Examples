test.scoregpu = function(){
  n = 20
  K = 10
  otime = rep(1,n)/seq(1,n)
  event = rep(c(0,1), n/2)
  expr = matrix(seq(1,K*n)/seq(K*n,1), K, n)

  sumstat = sum(sapply(1:K, function(k){summary(coxph(Surv(otime,event)~x, data=data.frame(x=rank(expr[k,]),otime, event)))[["sctest"]][1]}))
  checkEquals(sum(scoregpu(otime, event, expr, test="npcox", B = 0, index = FALSE, pval=FALSE)), sumstat, tolerance=0.001)
}
