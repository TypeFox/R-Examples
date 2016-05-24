### R code from vignette source 'cherry.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: NAEP
###################################################
library(cherry)
data(NAEP)


###################################################
### code chunk number 3: pickFisher
###################################################
pickFisher(NAEP, c("HI","MN","IA"))


###################################################
### code chunk number 4: pickFisher_all
###################################################
pickFisher(NAEP)


###################################################
### code chunk number 5: pickFisher_return
###################################################
res <- pickFisher(NAEP, silent=TRUE)
res


###################################################
### code chunk number 6: curveFisher
###################################################
res <- curveFisher(NAEP)
res


###################################################
### code chunk number 7: curveFisher_fig
###################################################
curveFisher(NAEP)


###################################################
### code chunk number 8: curveFisher_order
###################################################
curveFisher(NAEP, select=c(8,3,4,2), plot=FALSE)
curveFisher(NAEP, order=c(8,3,4,2), plot=FALSE)


###################################################
### code chunk number 9: pickSimes
###################################################
hom <- hommelFast(NAEP)
pickSimes(hom, c("HI","MN","IA"))
curveSimes(hom, plot=FALSE)


###################################################
### code chunk number 10: hommel
###################################################
hom <- hommelFast(NAEP, simes=FALSE)
pickSimes(hom, c("HI","MN","IA"))
curveSimes(hom, plot=FALSE)


###################################################
### code chunk number 11: hsu
###################################################
ps <- c(A = 0.051, B = 0.064, C = 0.097, D = 0.108)


###################################################
### code chunk number 12: fisher_local
###################################################
mytest <- function(hypotheses) {
  p.vals <- ps[hypotheses]
  m <- length(p.vals)
  statistic <- -2 * sum(log(p.vals))
  p.out <- pchisq(statistic, df=2*m, lower.tail=FALSE)
  return(p.out)
}


###################################################
### code chunk number 13: fisher_local_examples
###################################################
mytest("A")
mytest(c("B","C","D"))
mytest(names(ps))


###################################################
### code chunk number 14: Ftest
###################################################
hypotheses <- c("pop15", "pop75", "dpi", "ddpi")
fullfit <- lm(sr~., data=LifeCycleSavings)
myFtest <- function(hyps) {
  others <- setdiff(hypotheses, hyps)
  form <- formula(paste(c("sr~",  paste(c("1", others), collapse="+"))))
  anov <- anova(lm(form, data=LifeCycleSavings), fullfit, test="F")
  pvalue <- anov$"Pr("[2]        # NB replace Pr by P for for R < 2.14.0
  return(pvalue)
}
myFtest(c("pop15", "pop75"))


###################################################
### code chunk number 15: closed
###################################################
ct <- closed(mytest, names(ps))


###################################################
### code chunk number 16: closed_show
###################################################
ct


###################################################
### code chunk number 17: pick
###################################################
pick(ct, c("A", "B"))
pick(ct, c("C", "D"))


###################################################
### code chunk number 18: hypotheses
###################################################
hypotheses(ct)


###################################################
### code chunk number 19: pick_output
###################################################
res <- pick(ct, c("C", "D"), silent=TRUE)
res


###################################################
### code chunk number 20: defining
###################################################
defining(ct)


###################################################
### code chunk number 21: shortlist
###################################################
shortlist(ct)


###################################################
### code chunk number 22: closed_adjusted
###################################################
cta <- closed(mytest, names(ps), alpha = NA)


###################################################
### code chunk number 23: closed_adjusted2
###################################################
ctb <- closed(mytest, names(ps), alpha = 0.1, adjust = TRUE)


###################################################
### code chunk number 24: pick_adjusted
###################################################
pick(cta)


###################################################
### code chunk number 25: adjusted_function
###################################################
adjusted(cta, c("A", "B", "C"))
adjusted(cta, c("A", "B", "C"), n=1)


###################################################
### code chunk number 26: confidence_distribution
###################################################
pick(cta, names(ps), plot=TRUE)


###################################################
### code chunk number 27: defining_adjusted
###################################################
shortlist(cta, alpha=0.05)


###################################################
### code chunk number 28: alpha
###################################################
alpha(cta) <- 0.05


###################################################
### code chunk number 29: alphaNA
###################################################
alpha(cta) <- NA


###################################################
### code chunk number 30: generatedata
###################################################
set.seed(1)
n=100
p=20
X <- matrix(rnorm(n*p),n,p)
beta <- c(rep(0,2),rep(1,4),rep(0,2),rep(1,4),rep(0,2),
          rep(1,4),rep(0,2))
Y <- X %*% beta + rnorm(n)


###################################################
### code chunk number 31: Ftest
###################################################
mytest <- function(left,right)
{
  X <- X[,(left:right),drop=FALSE]
  lm.out <- lm(Y ~ X)
  x <- summary(lm.out)
  return(pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],
            lower.tail=FALSE))  
}


###################################################
### code chunk number 32: regionmethod
###################################################
reg <- regionmethod(rep(1,p), mytest, all_pvalues=TRUE, isadjusted=TRUE)
summary(reg)


###################################################
### code chunk number 33: implications
###################################################
implications(reg)


###################################################
### code chunk number 34: implications2
###################################################
implications(reg, alpha=0.01)


###################################################
### code chunk number 35: pvaluefunction
###################################################
pvalue(reg,1,20)


###################################################
### code chunk number 36: regionplot
###################################################
regionplot(reg)


###################################################
### code chunk number 37: regionplot2
###################################################
regionplot2(reg)


###################################################
### code chunk number 38: regionpick
###################################################
regionpick(reg, list(c(1,5),c(16,20)))


###################################################
### code chunk number 39: regionpick
###################################################
regionpick(reg, list(c(1,20)))


###################################################
### code chunk number 40: generatedataDAG
###################################################
set.seed(1)
n=100
p=4
X <- matrix(rnorm(n*p),n,p)
beta <- c(0,0.5,0.5,0)
Y <- X %*% beta + rnorm(n)


###################################################
### code chunk number 41: makesets
###################################################
sets <- list(c(1,2,3,4), c(1,2), c(2,3,4), c(2,3), 1, 2, 3, 4)
names(sets) <- c(1234, 12, 234, 23, 1, 2, 3, 4)


###################################################
### code chunk number 42: makeDAG
###################################################
struct <- construct(sets)


###################################################
### code chunk number 43: twoway
###################################################
istwoway(struct)


###################################################
### code chunk number 44: Ftest2
###################################################
mytest <- function(set)
{ 
  X <- X[,set,drop=FALSE]
  lm.out <- lm(Y ~ X)
  x <- summary(lm.out)
  return(pf(x$fstatistic[1],x$fstatistic[2],
            x$fstatistic[3],lower.tail=FALSE))  
}


###################################################
### code chunk number 45: DAGmethod
###################################################
DAG <- DAGmethod(struct, mytest, isadjusted=TRUE)
summary(DAG)


###################################################
### code chunk number 46: implicationsDAG
###################################################
implications(DAG)


###################################################
### code chunk number 47: pvaluefunctionDAG
###################################################
pvalue(DAG,4)


###################################################
### code chunk number 48: pvaluefunctionDAG2
###################################################
pvalue(DAG, "23")


###################################################
### code chunk number 49: DAGpick
###################################################
DAGpick(DAG, 1:3)


###################################################
### code chunk number 50: generatedataDAGholm
###################################################
set.seed(1)
n=100
p=4
X <- matrix(rnorm(n*p),n,p)
beta <- c(0,0.5,0.5,0)
Y <- X %*% beta + rnorm(n)


###################################################
### code chunk number 51: makesetsholm
###################################################
sets <- list(c(1,2,3,4), c(1,2), c(2,3,4), c(2,3), 1, 2, 3, 4)
names(sets) <- c(1234, 12, 234, 23, 1, 2, 3, 4)


###################################################
### code chunk number 52: makeDAGholm
###################################################
struct <- construct(sets)


###################################################
### code chunk number 53: twowayholm
###################################################
istwoway(struct)


###################################################
### code chunk number 54: Ftest2holm
###################################################
mytest <- function(set)
{ 
  X <- X[,set,drop=FALSE]
  lm.out <- lm(Y ~ X)
  x <- summary(lm.out)
  return(pf(x$fstatistic[1],x$fstatistic[2],
            x$fstatistic[3],lower.tail=FALSE))  
}


###################################################
### code chunk number 55: holmmethod
###################################################
DAG <- structuredHolm(struct, mytest, isadjusted=TRUE)
summary(DAG)


###################################################
### code chunk number 56: holmpval
###################################################
pvalue(DAG)


