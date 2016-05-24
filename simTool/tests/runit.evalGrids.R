test.evalGrids = function(){
  require(reshape)
  set.seed(15032013)
  dg = expandGrid(fun=c("runif"), n=1:4)  
  eg = evalGrids(dg)
  checkEquals(all(eg$dataGrid == dg), TRUE)
  checkEquals(all(eg$procGrid == expandGrid(proc="length")), TRUE)
  
  data = unlist(sapply(eg$simulation, function(l) l[[1]]$data))
  set.seed(15032013)
  checkEquals(!is.null(data), TRUE)
  checkEquals(all(data == runif(sum(1:4))), TRUE)
  
  results = unlist(sapply(eg$simulation, function(l) l[[1]]$results))
  checkEquals(all(results == 1:4), TRUE)
  
  set.seed(15032013)
  pg = expandGrid(proc=c("mean"))
  eg = evalGrids(dg, pg)

  lpkgN = function(dummy, name) c(ret=any(loadedNamespaces() == name))  
  pg = expandGrid(proc="lpkgN", name=c("boot", "MASS"))
  eg = evalGrids(dg, pg, ncpus=2, clusterLibraries=c("MASS", "boot"), rep=2, envir=environment())    
  checkEquals(all(as.data.frame(eg)$ret == TRUE), TRUE)
  
  ## summary.fun tests.  
  set.seed(23112013)
  dg = expandGrid(proc="runif", n=c(5, 100, 1000))
  pg = expandGrid(fun=c("summary", "mean"))
  checkException(evalGrids(dg, pg, post.proc=mean))  
  eg = evalGrids(dg, pg, replications=10)
  df1 = as.data.frame(eg, summary.fun=c(mean, sd))
  set.seed(23112013)
  eg = evalGrids(dg, pg, replications=10, summary.fun=c(mean, sd))
  df2 = as.data.frame(eg)
  df2$replication=NULL
  checkEquals(identical(df1, df2), TRUE)

  set.seed(23112013)
  dg = expandGrid(proc="runif", n=c(5, 100, 1000))
  pg = expandGrid(fun=c("summary", "mean"))
  eg = evalGrids(dg, pg, replications=10)
  postVec = function(results) c(mean=mean(results), sd=sd(results))
  df3 = as.data.frame(eg, summary.fun=postVec)
  checkEquals(identical(df3, df2), TRUE)
  set.seed(23112013)
  eg = evalGrids(dg, pg, replications=10, summary.fun=postVec)
  df2 = as.data.frame(eg)
  df2$replication=NULL
  checkEquals(identical(df3, df2), TRUE)
  

}
