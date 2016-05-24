test.as.data.frame.evalGrid = function(){
  dg = expandGrid(fun=c("runif"), n=1:4)
  pg = expandGrid(proc=c("mean"))  
  eg = evalGrids(dg, pg)
  df = as.data.frame(eg)
  tmp = c(
    eg$simulation[[1]][[1]]$results[[1]],
    eg$simulation[[2]][[1]]$results[[1]],
    eg$simulation[[3]][[1]]$results[[1]],
    eg$simulation[[4]][[1]]$results[[1]])
  
  checkEquals(all(df$V1==tmp), TRUE)

  checkException(as.data.frame(eg, value.fun=identity))
  checkException(as.data.frame(eg, post.proc=mean))
  
  # check if summary.fun can handle numeric and logical
  f = function(x) x <= 0.05
  df = as.data.frame(eg, convert.result.fun=f, summary.fun=c(mean, median))
  
  postVec = function(results) summary(results)
  df = as.data.frame(eg, convert.result.fun=f, summary.fun=postVec)
  
  f = function(x) c(length(x), min = min(x), max(x))
  pg = expandGrid(proc=c("f"))
  # must set environment, otherwise evalGrids will not find the function f()
  eg = evalGrids(dg, pg, envir=environment())
  df = as.data.frame(eg)

  tmp = rbind(
    eg$simulation[[1]][[1]]$results[[1]],
    eg$simulation[[2]][[1]]$results[[1]],
    eg$simulation[[3]][[1]]$results[[1]],
    eg$simulation[[4]][[1]]$results[[1]])
  
  checkEquals(all(names(df) == c("i", "j", "fun", "n", "proc", "replication", "V1", "min", "V2")), TRUE)  
  checkEquals(all(df[, c("V1", "min", "V2")]==tmp), TRUE)

  genRegData <- function(){
    data.frame(
      x = 1:10,
      y = rnorm(10, mean=1:10))
  }
  
  
  set.seed(19032013)
  eg <- evalGrids(
    expandGrid(fun="genRegData"),
    expandGrid(proc="lm", formula=c("y ~ x", "y ~ x + I(x^2)")),
    replications=10, envir=environment())

  lm2df = function(lm.object) {
    ret = coef(summary(lm.object))
    data.frame(covariable = rownames(ret), ret, check.names=FALSE)
  }
  df<-as.data.frame(eg, convert.result.fun=lm2df)  
  require(plyr)
  tmp = ldply(eg$simulation[[1]], function(v) rbind(lm2df(v$results[[1]]), lm2df(v$results[[2]])))
  tmp = arrange(tmp, Estimate)
  df = arrange(df, Estimate)
  checkEquals(all(df[,-(1:6)] == tmp), TRUE)

  
  df = as.data.frame(eg, convert.result.fun=lm2df, summary.fun=mean)
  require(reshape)
  tmp = ldply(eg$simulation[[1]], function(v) rbind(lm2df(v$results[[1]])))
  mtmp = melt(tmp, id=1)
  tmp1 = cast(mtmp, ... ~ variable, mean)
  
  tmp = ldply(eg$simulation[[1]], function(v) rbind(lm2df(v$results[[2]])))
  mtmp = melt(tmp, id=1)
  tmp2 = cast(mtmp, ... ~ variable, mean)
  
  tmp = rbind(tmp1, tmp2)
  checkEquals(all(df[,-(1:5)] == tmp), TRUE)
  
  
  df = as.data.frame(eg, convert.result.fun=lm2df, summary.fun=c(mean, sd))
  
  tmp = ldply(eg$simulation[[1]], function(v) rbind(lm2df(v$results[[1]])))
  mtmp = melt(tmp, id=1)
  tmp1 = cast(mtmp, ... ~ variable, c(mean, sd))

  tmp = ldply(eg$simulation[[1]], function(v) rbind(lm2df(v$results[[2]])))
  mtmp = melt(tmp, id=1)
  tmp2 = cast(mtmp, ... ~ variable, c(mean, sd))
  
  tmp = rbind(tmp1, tmp2)
  checkEquals(all(df[,-(1:5)] == tmp), TRUE)  
}

test.as.data.frame.from.fallback = function(){
  require(plyr)
  genData = function(N){N}
  brichtAb = function(data){    
    if (data == 3)
      stop("Bewusster Error: fallBackTest")
    1
  }
  
  dg = expandGrid(fun=c("genData"), N=1:5)  
  pg = expandGrid(proc=c("brichtAb"))  
  if(is.element("RUnitFallBack.Rdata", dir()))
    file.remove("RUnitFallBack.Rdata")
  
  checkException(eg <- evalGrids(dg, pg, replications=4, progress=TRUE, fallback="RUnitFallBack", envir=environment()))
  rm(list=ls())
  load("RUnitFallBack.Rdata")
  df = as.data.frame(fallBackObj)
  
  # dput
  df2 = structure(list(
    i = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 4L, 5L), 
    j = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), 
    fun = c("genData", "genData", "genData", "genData", 
            "genData", "genData", "genData", "genData", "genData", "genData", "genData"), 
    N = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 4L, 5L), 
    proc = c("brichtAb", "brichtAb", "brichtAb", "brichtAb", 
             "brichtAb", "brichtAb", "brichtAb", "brichtAb", 
             "brichtAb", "brichtAb", "brichtAb"), 
    replication = structure(c(1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, NA, NA, NA), 
                            .Label = c("1", "2", "3", "4"), class = "factor"), 
    V1 = c(1, 1, 1, 1, 1, 1, 1, 1, NA, NA, NA), 
    .evalGridComment = structure(c(NA, NA, NA, NA, NA, NA, NA, NA, 1L, 1L, 1L), 
                                 .Label = "Results missing", class = "factor")),
    .Names = c("i", "j", "fun", "N", "proc", "replication", "V1", ".evalGridComment"), 
    row.names = c(NA, -11L), class = "data.frame")
  checkEquals(df, df2)
}
