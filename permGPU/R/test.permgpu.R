test.permgpu = function(test){
  n = 20
  K = 10
  otime = rep(1,n)/seq(1,n)
  event = rep(c(0,1), n/2)
  pdat = data.frame(otime, event)
  expr = matrix(seq(1,K*n)/seq(K*n,1), K, n)

  rownames(expr) = paste("g", 1:K, sep="")
  patid = paste("id", 1:n, sep="")
  rownames(pdat) = patid
  colnames(expr) = patid

  testdat = makeExprSet(expr, pdat)

  if (test == "ttest"){
    sumstat = sum(sapply(1:K, function(k){t.test(x~y, data=data.frame(x=exprs(testdat)[k,], y=pData(testdat)$event))$statistic}))
    checkEquals(sum(permgpu(testdat, "event", B=0, test="ttest")[["stat"]]), sumstat, tolerance=0.001)
  }
  else if (test == "wilcoxon"){
    n0 = sum(event==0)
    n1 = sum(event==1)
    m  = n1 * (n0+n1+1.0) / 2.0
    std = sqrt(n0*n1*(n0+n1+1.0)/12.0)
    sumstat = sum(sapply(1:K, function(k){(sum(rank(exprs(testdat)[k,])[event==1])-m)/std}))
    checkEquals(sum(permgpu(testdat, "event", B=0, test="wilcoxon")[["stat"]]), sumstat, tolerance=0.001)
  }
  else if (test == "npcox"){
    sumstat = sum(sapply(1:K, function(k){summary(coxph(Surv(otime,event)~x, data=data.frame(x=rank(exprs(testdat)[k,]),otime=pData(testdat)[["otime"]], event=pData(testdat)[["event"]])))[["sctest"]][1]}))
    checkEquals(sum(permgpu(testdat, "otime", "event", B=0, test="npcox")[["stat"]]), sumstat, tolerance=0.001)
  }
  else if (test == "cox"){
    sumstat = sum(sapply(1:K, function(k){summary(coxph(Surv(otime,event)~x, data=data.frame(x=exprs(testdat)[k,],otime=pData(testdat)[["otime"]], event=pData(testdat)[["event"]])))[["sctest"]][1]}))
    checkEquals(sum(permgpu(testdat, "otime", "event", B=0, test="cox")[["stat"]]), sumstat, tolerance=0.001)
  }
  else
    stop("test unknown.")
}
