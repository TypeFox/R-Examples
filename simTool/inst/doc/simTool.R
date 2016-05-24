## ------------------------------------------------------------------------
library(simTool)
library(plyr)
library(reshape)
print(dg <- rbind.fill(
  expandGrid(fun="rexp", n=c(10, 20), rate=1:2),
  expandGrid(fun="rnorm", n=c(10, 20), mean=1:2)))

## ------------------------------------------------------------------------
print(pg<-rbind.fill(
  expandGrid(proc="min"),
  expandGrid(proc="mean", trim=c(0.1, 0.2))))

## ----, eval=FALSE--------------------------------------------------------
#  1.  convert dg to R-functions  {g_1, ..., g_k}
#  2.  convert pg to R-functions  {f_1, ..., f_L}
#  3.  initialize result object
#  4.  append dg and pg to the result object
#  5.  t1 = current.time()
#  6.  for g in  {g_1, ..., g_k}
#  7.      for r in 1:replications (optionally in a parallel manner)
#  8.          data = g()
#  9.          for f in  {f_1, \ldots, f_L}
#  10.             append f(data) to the result object
#  11.         optionally append data to the result object
#  12.      optionally summarize the result object over all
#           replications but separately for f_1, ..., f_L
#  13.     optionally save the results so far obtained to HDD
#  14. t2 = current.time()
#  15. Estimate the number of replications per hour from t1 and t2

## ------------------------------------------------------------------------
dg = expandGrid(fun="rnorm", n=10, mean=1:2)
pg = expandGrid(proc="min")
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 2)
as.data.frame(eg)

## ------------------------------------------------------------------------
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 3)
as.data.frame(eg)

## ------------------------------------------------------------------------
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 1000)
object.size(eg)
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 1000,
    discardGeneratedData = TRUE)
object.size(eg)

## ------------------------------------------------------------------------
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 10,
    progress = TRUE)

## ------------------------------------------------------------------------
dg = expandGrid(fun="runif", n=c(10,20,30))
pg = expandGrid(proc=c("min", "max"))
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 1000, 
    summary.fun=mean)
as.data.frame(eg)
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 1000, 
    summary.fun=c(mean, sd))
as.data.frame(eg)

## ------------------------------------------------------------------------
set.seed(1234)
# summarize the result objects as soon as possible
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 1000, 
    summary.fun=mean)
as.data.frame(eg)
set.seed(1234)
# keeping the result objects
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 1000)
# summarize the result objects by as.data.frame
as.data.frame(eg, summary.fun=mean)

## ------------------------------------------------------------------------
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 10, 
    ncpus=2, summary.fun=mean)
as.data.frame(eg)

## ------------------------------------------------------------------------
require(parallel)
cl = makeCluster(rep("localhost", 2), type="PSOCK") 
eg = evalGrids(dataGrid = dg, procGrid = pg, replications = 10, 
    cluster=cl, summary.fun=mean)
as.data.frame(eg)
stopCluster(cl)

## ------------------------------------------------------------------------
library(boot)
ratio <- function(d, w) sum(d$x * w)/sum(d$u * w)
city.boot <- boot(city, ratio, R = 999, stype = "w", 
    sim = "ordinary")
boot.ci(city.boot, conf = c(0.90, 0.95),
    type = c("norm", "basic", "perc", "bca"))

## ------------------------------------------------------------------------
returnCity = function(){
  city
}
bootConfInt = function(data){
city.boot <- boot(data, ratio, R = 999, stype = "w", 
    sim = "ordinary")
boot.ci(city.boot, conf = c(0.90, 0.95),
    type = c("norm", "basic", "perc", "bca"))  
}

## ------------------------------------------------------------------------
dg = expandGrid(fun="returnCity")
pg = expandGrid(proc="bootConfInt")
eg = evalGrids(dg, pg, replications=10, ncpus=2,
    clusterLibraries=c("boot"), 
    clusterGlobalObjects=c("ratio"))

## ------------------------------------------------------------------------
genData = function(n){
  n
}
anaData = function(data){
  if (data == 4) 
    stop("Simulated error that terminates the simulation")
  data^2
}
dg = expandGrid(fun="genData", n=1:5)
pg = expandGrid(proc="anaData")
try(eg <- evalGrids(dg, pg, replications=2, 
    fallback="simTool_fbTest"))

## ------------------------------------------------------------------------
# clean the current R-session
rm(list=ls())
load("simTool_fbTest.Rdata")
as.data.frame(fallBackObj)

## ------------------------------------------------------------------------
# masking summary from the base package
summary = function(x) sd(x) 
g = function(x) quantile(x, 0.1)
someFunc = function(){
  summary = function(x) c(sd=sd(x), mean=mean(x))
  
  dg = expandGrid(fun="runif", n=100)
  pg = expandGrid(proc=c("summary", "g"))

  # the standard is to use the global
  # environment, hence summary defined outside
  # of someFunc() will be used
  print(as.data.frame(evalGrids(dg, pg)))
  cat("--------------------------------------------------\n")
  # will use the local defined summary, but g
  # from the global environment, because
  # g is not locally defined.
  print(as.data.frame(evalGrids(dg, pg, envir=environment())))
}
someFunc()

## ------------------------------------------------------------------------
dg = rbind.fill(
  expandGrid(fun="rexp", n=c(10, 20), rate=1:2),
  expandGrid(fun="rnorm", n=c(10, 20), mean=1:2))
pg = rbind.fill(
  expandGrid(proc="min"),
  expandGrid(proc="mean", trim=c(0.1, 0.2)))

## ------------------------------------------------------------------------
eg = evalGrids(dg, pg, replications=100)

## ------------------------------------------------------------------------
names(eg)

## ------------------------------------------------------------------------
dg[7,]

## ------------------------------------------------------------------------
pg[3,]

## ------------------------------------------------------------------------
eg$simulation[[7]][[22]]$results[[3]]

## ------------------------------------------------------------------------
mean(eg$simulation[[7]][[22]]$data, trim=0.2)

## ------------------------------------------------------------------------
genRegData <- function(){
  data.frame(
    x = 1:10,
    y = rnorm(10, mean=1:10))
}

## ------------------------------------------------------------------------
eg <- evalGrids(
  expandGrid(fun="genRegData"),
  expandGrid(proc="lm", formula=c("y ~ x", "y ~ x + I(x^2)")),
  replications=100)
class(eg$simulation[[1]][[1]]$results[[1]])

## ------------------------------------------------------------------------
head(df<-as.data.frame(eg, convert.result.fun=coef))

## ------------------------------------------------------------------------
as.data.frame(eg, convert.result.fun=coef, summary.fun=c(mean, sd))

