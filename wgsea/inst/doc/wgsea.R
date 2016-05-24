### R code from vignette source 'wgsea.Rnw'

###################################################
### code chunk number 1: wgsea.Rnw:64-78
###################################################
library(snpStats)
library(wgsea)

## load example data from snpStats
data(for.exercise,package="snpStats")

## generate an artificial indicator of test and control SNPs.
snpsum <- col.summary(snps.10)
snp.indicator <- abs(snpsum$z.HWE) > 1.96

## subset to case and control objects
case <- snps.10[subject.support$cc==1,!is.na(snp.indicator)]
control <- snps.10[subject.support$cc==0,!is.na(snp.indicator)]
snp.indicator <- snp.indicator[!is.na(snp.indicator)]


###################################################
### code chunk number 2: wgsea.Rnw:86-88
###################################################
maf <- col.summary(control)[,"MAF"]
p <- pairtest(case,control)


###################################################
### code chunk number 3: wgsea.Rnw:92-94
###################################################
wilcox.test(x=p[snp.indicator==1], 
            y=p[snp.indicator==0])


###################################################
### code chunk number 4: wgsea.Rnw:109-111
###################################################
n.perm=10
p.perm <- pairtest(case,control,n.perm=n.perm)


###################################################
### code chunk number 5: wgsea.Rnw:117-120
###################################################
## calculate wilcoxon score
W <- wilcoxon(p,snps.in=which(snp.indicator==1))
Wstar <- wilcoxon(p.perm,snps.in=which(snp.indicator==1))


###################################################
### code chunk number 6: wgsea.Rnw:143-157
###################################################
## running empirical estimate of variance
#W.var <- sapply(2:n.perm,function(i) var(Wstar[sample(1:n.perm,size=i)]))
W.var <- sapply(2:n.perm,function(i) var(Wstar[1:i]))

## theoretical variance for comparison
n1 <- sum(snp.indicator)
n2 <- sum(!snp.indicator)
var.theoretical <-  exp(log(n1) + log(n2) + log(n1+n2+1) - log(12))

## plot
plot(2:n.perm,W.var,ylim=range(c(W.var,var.theoretical)),pch=".",
     xlab="Permutation number",ylab="Var(W*)",main="Estimate of Var(W*) vs number of permutations",
     sub="(dotted line shows theoretical value)")
abline(h=var.theoretical,lty=3)


###################################################
### code chunk number 7: wgsea.Rnw:178-179
###################################################
Z.value(W=W, Wstar=Wstar, n.in=sum(snp.indicator==1), n.out=sum(snp.indicator==0))


###################################################
### code chunk number 8: wgsea.Rnw:190-191
###################################################
cor.test(maf,p)


###################################################
### code chunk number 9: wgsea.Rnw:205-210
###################################################
## calculate wilcoxon score
W.weighted <- wilcoxon(p,snps.in=which(snp.indicator==1),weights=maf)
Wstar.weighted <- wilcoxon(p.perm,snps.in=which(snp.indicator==1),weights=maf)
Z.value(W=W.weighted, Wstar=Wstar.weighted,
        n.in=sum(snp.indicator==1), n.out=sum(snp.indicator==0))


###################################################
### code chunk number 10: wgsea.Rnw:247-250
###################################################
Z.value(W=list(W,W), Wstar=list(Wstar,Wstar), 
        n.in=rep(sum(snp.indicator==1),2), 
        n.out=rep(sum(snp.indicator==0),2))


