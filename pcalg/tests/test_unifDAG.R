library(pcalg)
## setwd("/sfs/u/kalischm/research/packages/unifDAGs/")
## source("aux_general.R")
## source("unifDAG.R")
## source("unifDAG_approx.R")

## exact
set.seed(123)
n <- 2
reps <- 100
res <- matrix(0, reps, 3)
for (i in 1:reps) {
    gAM <- unifDAG(n)
    wm <- wgtMatrix(gAM)
    if (wm[1,2] != 0) {
        res[i,1] <- 1
    } else {
          if (wm[2,1] != 0) {
              res[i,2] <- 1
          } else {
                res[i,3] <- 1
            }
      }
}
if (!all(round(colMeans(res), 2) == c(0.41, 0.31, 0.28))) {
    stop("test_unifDAG: Saved result could not be reproduced (unifDAG).")
}
## apply(res,2,function(x) sd(x)/sqrt(length(x)))
## Result for reps = 10'000
## 0.3409 0.3284 0.3307
## 0.004740355 0.004696547 0.004704887

## Approx
set.seed(123)
n <- 2
reps <- 100
res <- matrix(0, reps, 3)
for (i in 1:reps) {
    gAM <- unifDAG.approx(n = 2, n.exact = 2)
    wm <- wgtMatrix(gAM)
    if (wm[1,2] != 0) {
        res[i,1] <- 1
    } else {
          if (wm[2,1] != 0) {
              res[i,2] <- 1
          } else {
                res[i,3] <- 1
            }
      }
}
colMeans(res)
## apply(res,2,function(x) sd(x)/sqrt(length(x)))
## Result for reps = 10'000
## 0.3332 0.3401 0.3267
## 0.004713809 0.004737662 0.004690300

if (!all(round(colMeans(res), 2) == c(0.34, 0.31, 0.35))) {
    stop("test_unifDAG: Saved result could not be reproduced (unifDAG.approx).")
}

## check weights
set.seed(123)
n <- 120
g <- unifDAG.approx(n=n, n.exact = 15, weighted = TRUE, wFUN = list(runif, min=0, max=1))
m <- wgtMatrix(g)
v <- as.numeric(m)
v <- v[v!=0]
ct <- cut(x=v, breaks=seq(0,1,by=0.1))
pval <- chisq.test(as.numeric(table(ct)), p = rep(0.1,10))$p.value
if (pval < 0.05) {
    stop("test_unifDAG: Edge weights don't look uniformly distributed!")
}

