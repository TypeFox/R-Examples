library(pcalg)

## Create a weighted DAG --- Example from ../man/jointIda.Rd  (keep in sync !)
p <- 6
V <- as.character(1:p)
edL <- list(
  "1" = list(edges=c(3,4), weights=c(1.1,0.3)),
  "2" = list(edges=c(6),  weights=c(0.4)),
  "3" = list(edges=c(2,4,6),weights=c(0.6,0.8,0.9)),
  "4" = list(edges=c(2),weights=c(0.5)),
  "5" = list(edges=c(1,4),weights=c(0.2,0.7)),
  "6" = NULL)
myDAG <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed") ## true DAG
myCPDAG <- dag2cpdag(myDAG) ## true CPDAG
covTrue <- trueCov(myDAG) ## true covariance matrix

n <- 1000
## simulate Gaussian data from the true DAG
dat <- if (require("mvtnorm")) {
  set.seed(123)
  rmvnorm(n, mean=rep(0,p), sigma=covTrue)
} else readRDS(system.file(package="pcalg", "external", "N_6_1000.rds"))

## estimate CPDAG -- see  help(pc)
suffStat <- list(C = cor(dat), n = n)
pc.fit <- pc(suffStat, indepTest = gaussCItest, p = p, alpha = 0.01, u2pd="relaxed")

## Suppose that we know the true CPDAG and covariance matrix
m1 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myCPDAG, technique="RRC")
m2 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, graphEst=myCPDAG, technique="MCD")

## Instead of knowing the true CPDAG, it is enough to know only
## all jointly valid parent sets of the intervention variables ## to use RRC or MCD
ajv.pasets <- list(list(5, c(3,4)),
                   list(integer(0), c(3,4)),
                   list(3, c(3,4)))
m3 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, all.pasets=ajv.pasets, technique="RRC")
m4 <- jointIda(x.pos=c(1,2), y.pos=6, covTrue, all.pasets=ajv.pasets, technique="MCD")

## MM: FIXME all(. == .) for  0.99 !!!!  ---- and similarly below !
mTrue <- matrix(c(0.99, 0.4, 0.99, 0.4, 0, 0.4), 2,3)
Rnd <- 7
res1 <- (all(round(m1,Rnd) == mTrue) &
         all(round(m2,Rnd) == mTrue) &
         all(round(m3,Rnd) == mTrue) &
         all(round(m4,Rnd) == mTrue) )

if(!res1) stop("Test in jointIda: True causal effects were not recovered!")

## From the true DAG, we can compute the true total joint effects
## using RRC or MCD
jointIda(x.pos=c(1,2),y.pos=6,covTrue,graphEst=myDAG,technique="RRC")
jointIda(x.pos=c(1,2),y.pos=6,covTrue,graphEst=myDAG,technique="MCD")

## jointIda also works when x.pos has length 1 and in the following example
## it gives the same result as ida() (see Note)
##
## When the CPDAG is known
v1 <- round(jointIda(x.pos=1,y.pos=6,covTrue,graphEst=myCPDAG,
                     technique="RRC"), Rnd)
v2 <- round(jointIda(x.pos=1,y.pos=6,covTrue,graphEst=myCPDAG,
               technique="MCD"), Rnd)
v3 <- round(ida(x.pos=1,y.pos=6,covTrue,graphEst=myCPDAG,
          method="global"), Rnd)
res2 <- (all(v1==v2) & all(v2==v3) & all(v3==v1))

if(!res2) stop("Test in jointIda: ida() and jointIda() don't produce the same result!")


