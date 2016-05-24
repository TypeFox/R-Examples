## ----setup, include=FALSE------------------------------------------------

knitr::opts_chunk$set(fig.pos="h")
knitr::opts_chunk$set(cache.path='./CovarianceEstimation_cache/')

## ----load_packages, results='hide', echo=FALSE, warning=FALSE, error=FALSE, message=FALSE----
library(knitcitations)
cleanbib()
options("citation_format" = "pandoc")

library(covmat)
library(xts)
library(robust)
library(PortfolioAnalytics)
library(rmgarch)

## ----load, eval=FALSE----------------------------------------------------
#  library(devtools)
#  install_github("arorar/covmat")

## ----doi, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE----------
bib <- read.bibtex("references.bib")

## ----rmt-load------------------------------------------------------------

data("dow30data")


## ---- eval=FALSE---------------------------------------------------------
#  
#  estRMT(R, Q =NA, cutoff = c("max", "each"),
#         eigenTreat = c("average", "delete") ,
#         numEig=1, parallel = TRUE)
#  

## ----rmt-est, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE------
model <- estRMT(dow30data, parallel=FALSE)


## ----rmt-plot, fig.width=8, fig.height=4, fig.keep='last'----------------

plot(model)


## ----rmt-custommoment----------------------------------------------------

custom.portfolio.moments <- function(R, portfolio) {
  momentargs <- list()
  momentargs$mu  <-  matrix(as.vector(apply(R,2, "mean")), ncol = 1)
  momentargs$sigma  <-  estRMT(R, parallel=FALSE)$cov
  momentargs$m3 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^2)
  momentargs$m4 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^3)

  return(momentargs)
}

## ----rmt-portfoliospec---------------------------------------------------

pspec.lo <- portfolio.spec(assets = colnames(dow30data))

#long-only
pspec.lo <- add.constraint(pspec.lo, type="full_investment")
pspec.lo <- add.constraint(pspec.lo, type="long_only")

pspec.lo <- add.objective(portfolio=pspec.lo, type="return", name="mean")
pspec.lo <- add.objective(portfolio=pspec.lo, type="risk", name="var")


## ----rmt-run, warning=FALSE, error=FALSE, message=FALSE, eval = FALSE----
#  
#  opt.ordinary <-
#    optimize.portfolio.rebalancing(dow30data, pspec.lo,
#                                   optimize_method="quadprog",
#                                   rebalance_on="months",
#                                   training_period=120,
#                                   trailing_periods=120)
#  opt.rmt <-
#    optimize.portfolio.rebalancing(dow30data, pspec.lo,
#                                   optimize_method="quadprog",
#                                   momentFUN = "custom.portfolio.moments",
#                                   rebalance_on="months",
#                                   training_period=120,
#                                   trailing_periods=120)

## ----rmt-results, eval = FALSE-------------------------------------------
#  
#  ordinary.wts <- na.omit(extractWeights(opt.ordinary))
#  ordinary <- Return.rebalancing(R=dow30data, weights=ordinary.wts)
#  
#  rmt.wts <- na.omit(extractWeights(opt.rmt))
#  rmt <- Return.rebalancing(R=dow30data, weights=rmt.wts)
#  
#  rmt.strat.rets <- merge.zoo(ordinary,rmt)
#  colnames(rmt.strat.rets) <- c("ordinary", "rmt")
#  

## ----rmt-results-main,  results='hide', echo=FALSE, warning=FALSE, error=FALSE, message=FALSE----

filepath <- "./CovarianceEstimation_cache/rmt_strategy_rets.RData"

if (file.exists(filepath)) {
  load(filepath)
} else {
  opt.ordinary <- 
  optimize.portfolio.rebalancing(dow30data, pspec.lo, 
                                 optimize_method="quadprog",
                                 rebalance_on="months", 
                                 training_period=120,
                                 trailing_periods=120)
  opt.rmt <- 
    optimize.portfolio.rebalancing(dow30data, pspec.lo, 
                                   optimize_method="quadprog",
                                   momentFUN = "custom.portfolio.moments",
                                   rebalance_on="months", 
                                   training_period=120, 
                                   trailing_periods=120)
  
  ordinary.wts <- na.omit(extractWeights(opt.ordinary))
  ordinary <- Return.rebalancing(R=dow30data, weights=ordinary.wts)
  
  rmt.wts <- na.omit(extractWeights(opt.rmt))
  rmt <- Return.rebalancing(R=dow30data, weights=rmt.wts)
  
  rmt.strat.rets <- merge.zoo(ordinary,rmt)
  colnames(rmt.strat.rets) <- c("ordinary", "rmt")
  
  save(rmt.strat.rets, file = filepath)
}
                              

## ----rmtstratplots-------------------------------------------------------
charts.PerformanceSummary(rmt.strat.rets,wealth.index = T, 
                          colorset = c("red","blue"), 
                          main=paste(c("Comparison of Portflio ",
                                     "Performance using two ",
                                     "different covariance matrices"),
                                     collapse=""), cex.legend = 1.3, 
                          cex.axis = 1.3, legend.loc = "topleft")


## ----isdcc-load----------------------------------------------------------

data("etfdata")


## ---- eval=FALSE---------------------------------------------------------
#  
#  isdccfit(R, numRegimes, transMatbounds = c(2,10),
#           dccBounds = c(0,1), w = NA, ...)
#  

## ----isdcc-est, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE, eval=FALSE----
#  model.isdcc <- isdccfit(etfdata, numRegimes=3, parallelType = 1, itermax = 100)
#  

## ----isdcc-est-main, results='hide' ,echo=FALSE, warning=FALSE, error=FALSE, message=FALSE----

filepath <- "./CovarianceEstimation_cache/isdcc_model.RData"
if (file.exists(filepath)) {
  load(filepath)
} else {
  model.isdcc <- isdccfit(etfdata, numRegimes=3, parallelType = 1, itermax = 100)
  save(model.isdcc, file = filepath)
}


## ----isdcc-plot, fig.width=8, fig.height=4, fig.keep='last'--------------

plot(model.isdcc)


## ----isdcc-custommoment--------------------------------------------------

custom.portfolio.moments.isdcc <- function(R, portfolio) {
  
  momentargs <- list()
  momentargs$mu  <-  matrix(as.vector(apply(R,2, "mean")), ncol = 1)
  
  result <- isdccfit(R, numRegimes = 2, itermax=100)
  ldate <- as.character(tail(index(R),1))
  maxProbIndex <- which.max(result$filtProb[ldate,])
  momentargs$sigma  <-  result$cov[[ldate]][[maxProbIndex]]
  
  momentargs$m3 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^2)
  momentargs$m4 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^3)
  
  return(momentargs)
}



custom.portfolio.moments.dcc <- function(R, portfolio) {
  
  momentargs <- list()
  momentargs$mu  <-  matrix(as.vector(apply(R,2, "mean")), ncol = 1)
  
  garch11.spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)), 
                            variance.model = list(garchOrder = c(1,1)
                                          , model = "sGARCH"), 
                            distribution.model = "norm")
  
  dcc.garch11.spec <- dccspec(uspec = multispec( replicate(ncol(R), 
                                                  garch11.spec) ), 
                             dccOrder = c(1,1), distribution = "mvnorm")
  
  dcc.fit <- dccfit(dcc.garch11.spec, data = R)
  momentargs$sigma  <-  rcov(dcc.fit)[,,as.character(tail(index(R),1))]
  
  momentargs$m3 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^2)
  momentargs$m4 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^3)
  
  return(momentargs)
}


## ----isdcc-portfoliospec-------------------------------------------------

datap <- etfdata["2009-07/"]

pspec.lo.isdcc <- portfolio.spec(assets = colnames(datap))

#long-only
pspec.lo.isdcc <- add.constraint(pspec.lo.isdcc, type="full_investment")
pspec.lo.isdcc <- add.constraint(pspec.lo.isdcc, type="long_only")

pspec.lo.isdcc <- add.objective(portfolio=pspec.lo.isdcc, 
                                type="return", name="mean")
pspec.lo.isdcc <- add.objective(portfolio=pspec.lo.isdcc, 
                                type="risk", name="var")


## ----isdcc-run, warning=FALSE, error=FALSE , message=FALSE, results='hide', eval=FALSE----
#  
#  ordinary <-
#    optimize.portfolio.rebalancing(datap, pspec.lo.isdcc,
#                                   optimize_method="quadprog",
#                                   rebalance_on="months",
#                                   training_period=120,
#                                   trailing_periods=120)
#  
#  opt.dcc <-
#    optimize.portfolio.rebalancing(datap, pspec.lo.isdcc,
#                                    optimize_method="quadprog",
#                                    momentFUN =
#                                     "custom.portfolio.moments.dcc",
#                                    rebalance_on="months",
#                                    training_period=120,
#                                    trailing_periods=120)
#  
#  opt.isdcc <-
#    optimize.portfolio.rebalancing(datap, pspec.lo.isdcc,
#                                            optimize_method="quadprog",
#                                            momentFUN =
#                                     "custom.portfolio.moments.isdcc",
#                                            rebalance_on="months",
#                                            training_period=120,
#                                            trailing_periods=120)

## ----isdcc-results,  eval=FALSE------------------------------------------
#  
#  ord.wts <- na.omit(extractWeights(ordinary))
#  ord <- Return.rebalancing(R=datap, weights=ord.wts)
#  
#  dcc.wts <- na.omit(extractWeights(opt.dcc))
#  dcc <- Return.rebalancing(R=datap, weights=dcc.wts)
#  
#  isdcc.wts <- na.omit(extractWeights(opt.isdcc))
#  isdcc <- Return.rebalancing(R=datap, weights=isdcc.wts)
#  
#  isdcc.strat.rets <- merge.zoo(merge.zoo(ord, dcc), isdcc)
#  colnames(isdcc.strat.rets) <- c("ordinary", "dcc", "isdcc")
#  

## ----isdcc-results-main,  results='hide', echo=FALSE, warning=FALSE, error=FALSE, message=FALSE----

filepath <- "./CovarianceEstimation_cache/isdcc_strategy_rets.RData"
if (file.exists(filepath)) {
  load(filepath)
} else {
  ordinary <- 
  optimize.portfolio.rebalancing(datap, pspec.lo.isdcc, 
                                 optimize_method="quadprog",
                                 rebalance_on="months", 
                                 training_period=120,
                                 trailing_periods=120)

  opt.dcc <- 
    optimize.portfolio.rebalancing(datap, pspec.lo.isdcc, 
                                    optimize_method="quadprog",
                                    momentFUN = 
                                     "custom.portfolio.moments.dcc",
                                    rebalance_on="months",
                                    training_period=120,
                                    trailing_periods=120)

  opt.isdcc <- 
    optimize.portfolio.rebalancing(datap, pspec.lo.isdcc, 
                                            optimize_method="quadprog",
                                            momentFUN = 
                                     "custom.portfolio.moments.isdcc",
                                            rebalance_on="months",
                                            training_period=120,
                                            trailing_periods=120)

  ord.wts <- na.omit(extractWeights(ordinary))
  ord <- Return.rebalancing(R=datap, weights=ord.wts)
  
  dcc.wts <- na.omit(extractWeights(opt.dcc))
  dcc <- Return.rebalancing(R=datap, weights=dcc.wts)
  
  isdcc.wts <- na.omit(extractWeights(opt.isdcc))
  isdcc <- Return.rebalancing(R=datap, weights=isdcc.wts)
  
  isdcc.strat.rets <- merge.zoo(merge.zoo(ord, dcc), isdcc)
  colnames(isdcc.strat.rets) <- c("ordinary", "dcc", "isdcc")
  
  save(isdcc.strat.rets, file = filepath)
}


## ----isdccstratplots, fig.width=10, fig.height=8-------------------------
charts.PerformanceSummary(isdcc.strat.rets,wealth.index = T,
                          colorset = c("red","blue","green"), 
                          main=paste(c("Comparison of Portflio ",
                                     "Performance using two ",
                                     "different covariance matrices"),
                                     collapse=""), cex.legend = 1.3, 
                          cex.axis = 1.3, legend.loc = "topright") 

## ----rmtdata-load--------------------------------------------------------

data("rmtdata")


## ---- eval=FALSE---------------------------------------------------------
#  
#  estSpikedCovariance(R, gamma = NA,
#                        numOfSpikes = NA,
#                        method = c("KNTest", "median-fitting"),
#                        norm = c("Frobenius", "Operator", "Nuclear"),
#                        pivot = 1, statistical = NA,
#                        fit = NA)
#  

## ----spike-plot, eval=FALSE----------------------------------------------
#  
#    plotSpikedCovariance(rmtdata)
#  

## ----spike-plot-main,  results='hide',  echo=FALSE, warning=FALSE, error=FALSE, message=FALSE, fig.show='hide'----

filepath <- "./CovarianceEstimation_cache/Shrinkage.png"
if (!file.exists(filepath)) {
  png(file = filepath, width = 12, height = 10, units = "in", res = 300)
  plotSpikedCovariance(rmtdata)
  dev.off()
}


## ---- eval=FALSE---------------------------------------------------------
#  
#  robustMultExpSmoothing(R, smoothMat = NA, startup_period = 10,
#                           training_period = 60 , seed = 9999, trials = 50,
#                           method = "L-BFGS-B", lambda = 0.2)
#  

## ----croux-est, eval=FALSE-----------------------------------------------
#  
#  data("dow30data")
#  symbols <- c('AAPL', 'GE', 'MSFT' , 'NKE', 'V')
#  
#  R <- dow30data[,which(colnames(dow30data) %in% symbols)]
#  smoothfit <- robustMultExpSmoothing(R)
#  

## ----croux-est-main, results='hide', echo=FALSE, warning=FALSE, error=FALSE, message=FALSE----
library(optimx)

data("dow30data")
symbols <- c('AAPL', 'GE', 'MSFT' , 'NKE', 'V')

R <- dow30data[,which(colnames(dow30data) %in% symbols)]

filepath <- "./CovarianceEstimation_cache/croux_smoothfit.RData"
if (file.exists(filepath)) {
  load(filepath)
} else {
  smoothfit <- robustMultExpSmoothing(R) 
  save(smoothfit, file = filepath)
}


## ----croux-plots1, fig.width=10, fig.height=6----------------------------

plotmissing(R)


## ----croux-plots2, fig.width=10, fig.height=6----------------------------

plotmissing(smoothfit$cleanValues)


## ----croux-cov-plots, fig.width=10, fig.height=6-------------------------

rob <- covRob(R)$cov
compareCov(smoothfit$covMat, rob, labels = c("Robust Croux", "Robust"))



## ----references, echo=FALSE, warning=FALSE, error=FALSE, message=FALSE, echo=FALSE----
#write.bibtex(file="references.bib")

