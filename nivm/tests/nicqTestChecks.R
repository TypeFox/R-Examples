library(nivm)
## following gives diff in prop outside CI.
## This is because ic!=ceiling(nc*q)
## Now gives warning
#x<-nicqTest(20,g=nimDiffOR,delta0=.1,q=.2,nc=200,nt=300,
#     ic=round(600*.2),conf.int=TRUE)
#x
## check that it works without specifying ic
#x<-nicqTest(20,g=nimDiffOR,delta0=.1,q=.2,nc=200,
#     nt=300,conf.int=TRUE)
#x
## check that alternative="greater" works
## x=114 barely rejects at 0.025 level
#x<-nicqTest(114,g=nimDiffOR,delta0=.1,q=.2,nc=200,
#      nt=300,conf.int=TRUE,alternative="greater")
## x=113 barely fails to reject at 0.025 level
#x<-nicqTest(113,g=nimDiffOR,delta0=.1,q=.2,nc=200,
#      nt=300,conf.int=TRUE,alternative="greater")