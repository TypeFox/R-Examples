## Test that perm gives the same answers in some special cases
library(perm)
library(coin)

#########################
## Wilcoxon Rank Sum Test
#########################
set.seed(1)
y<-rpois(10,5)
## Try and example with ties
table(y)
g<-as.factor(c(rep(1,4),rep(2,6)))

## Asymptotic method is the uncorrected one in wilcox.test
permTS(rank(y)~g,method="pclt")
wilcox.test(y~g,correct=FALSE)

## Compare to exact method in coin
permTS(rank(y)~g,method="exact.ce",alternative="two.sided")
permTS(rank(y)~g,method="exact.network",alternative="two.sided")

permTS(rank(y)~g,method="exact.ce",alternative="two.sided", control=permControl(tsmethod="abs"))
permTS(rank(y)~g,method="exact.network",alternative="two.sided",control=permControl(tsmethod="abs"))

## Note that coin uses the default two.sided p-value that matches our 
## alternative="two.sided" with permControl(tsmethod="abs"), but the default is permControl(tsmethod="central")
## need to use coin because wilcox.test exact=TRUE does not handle ties
wilcox_test(y~g,distribution=exact())

## Try one-sided
permTS(rank(y)~g,method="exact.network",alternative="less")
wilcox_test(y~g,distribution=exact(),alternative="less")

#################################
# Kruskal-Wallis Test
################################
g2<-as.factor(c(rep(1,3),rep(2,4),rep(3,3)))
kruskal.test(y~g2)
permKS(rank(y)~g2,exact=FALSE)

## Since both coin and perm use Monte Carlo, 
## we cannot expect an exact match
set.seed(11)
kruskal_test(y~g2,distribution="approximate")
permKS(rank(y),g2,method="exact.mc")

##############################
# Trend Tests
###############################
g3<-as.numeric(g2)
independence_test(y~g3,distribution="asymptotic",teststat="max")
independence_test(y~g3,distribution="asymptotic",teststat="quad")
independence_test(y~g3,distribution="asymptotic",teststat="scalar")
permTREND(y,g3,method="pclt")

independence_test(y~g3,alternative="less",distribution="asymptotic",teststat="max")
independence_test(y~g3,alternative="less",distribution="asymptotic",teststat="scalar")
permTREND(y,g3,alternative="less",method="pclt")

# I tested this data set using the Linear-by-Linear 
# Association test in StatXact Procs 8, the Asymptotic 
# p-value matches the two-sided one from permTREND
# asymptotic p=0.8604
# exact p = 0.9305
#permTREND(y,g3,method="exact.mc",alternative="less",nmc=10^5-1)
# gives p=.5945 with 99 pcnt ci (.5905,.5985)
#independence_test(y~g3,alternative="less",distribution=approximate(10^5-1),teststat="scalar")
# gives p=.593 
