## This is a hard two sample precision problem.
## What makes it hard is that the sum of the 
## first 2 values of y are the same as the sum
## of the second two, and to get the correct 
## answer you need to get them exactly equal
## You should get the same answer regardless of 
## whether you multiply y by a constant
## Here both coin 1.0-8 and perm 0.9.1 give the correct answer
## (earlier versions of coin had problems with this example)
library(coin)
packageDescription("coin")$Version
library(perm)
y<-c(1.11,2.22,3.33,0,100,1)
g<-as.factor(c(1,1,0,0,1,0))
independence_test(y~g,distribution=exact())
permTS(y~g,method="exact.ce",alternative="two.sided",control=permControl(tsmethod="abs"))
independence_test(y~g,distribution=exact(),alternative="less")
permTS(y~g,method="exact.ce",alternative="less")
independence_test(y~g,distribution=exact(),alternative="greater")
permTS(y~g,method="exact.ce",alternative="greater")

y<-4.44*y
independence_test(y~g,distribution=exact())
permTS(y~g,method="exact.ce",alternative="two.sided",control=permControl(tsmethod="abs"))