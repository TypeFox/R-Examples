######################################################################
## File: test1.R
## Author: Ray Griner
## Date: November 11, 2011 (111111!)
## Desc: Short examples of predicted interval plots (PIPs)
######################################################################

set.seed(12345)

#######################################################
# EXAMPLE 1: TWO SAMPLE, BINARY OUTCOME
#######################################################
myY<-c(rep(1,times=20),rep(0,times=80),rep(1,times=25),rep(0,times=25))
myGroup<-c(rep('A',100),rep('B',50))

## Predicted intervals near the 50th percentile are darkest, and get
##  lighter as they go towards the 0th and 100th percentile
my.col.fun<-function(dist) { gray( 6/5*abs(dist/100-.5)+.3) }

pips <- pred.int(y=myY, group=myGroup, N=c(400,400), data.type="binary", iters=100)

print(pips, pi.count=100)
plot(pips, pi.col.fun=my.col.fun)


