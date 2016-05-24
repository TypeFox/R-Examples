## test the different arguments in the nicq function
set.seed(1)
Nc<-100
Nt<-200
Yt<-rpois(Nt,10)
Yc<-rpois(Nc,10)
Y<-c(Yt,Yc)
Z<-c(rep(2,Nt),rep(1,Nc))

## Show that all three formats give the same answer: 

out1<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z)
out2<-nicqTest(x=Yt,g=nimDiffOR,delta0=.1,q=.2,yc=Yc)
out3<-nicqTest(x=out2$statistic,g=nimDiffOR,delta0=.1,q=.2,ic=out2$parameter["i"],nc=Nc,nt=Nt)


all.equal(out1,out2)
all.equal(out1,out3)
## method description different for ties, because input.type=3 assumes no ties
out1$method
out3$method

## make sure censoring after the ith control failure does not change the results
T1i<-sort(Yc)[out2$parameter["i"]]
Status<-rep(1,300)
Status[Y>T1i]<-0
out1b<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z,status=Status)
all.equal(out1,out1b)
## but censoring at T1i should give an error
#Status[Y>=T1i]<-0
#out1b2<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z,status=Status)



## check for p-values matching the confidence intervals 
## set conf.level = 1-2*p.value so that the upper limit 
## should equal the null difference (delta0) value
out1c<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z,conf.level=1-2*out1$p.value)
out1c$conf.int
out1c$null.value
#all.equal(out1,out1c)

## check ties="cons". F1(t)=0.20 is at t=7, and ties at yc=7, so be conservative 
## and take x2(T1i)=all the failed at or before yc=7
out1$parameter["i"]
cumsum(table(Yc))
cumsum(table(Yt))
out1$statistic

## for Poisson example, when ties="approx" then i=
##  x2T1i=round(33+ (20-17)*(51-33)/(24-17+1)) = 40
## with large numbers of ties, this makes a big difference in the p-value 
out4<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z,ties="approx")
out4$statistic
out4

#### Also can rerun with normal data to make sure continuous data work
set.seed(3)
Nc<-10
Nt<-20
Yt<-rnorm(Nt)
Yc<-rnorm(Nc)
Y<-c(Yt,Yc)
Z<-c(rep(2,Nt),rep(1,Nc))
## use small sample size to 
## see how much the g function matters: 

out5<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z)
out6<-nicqTest(x=Y,g=nimDiff,delta0=.1,q=.2,z=Z)
out7<-nicqTest(x=Y,g=nimOR,delta0=.1,q=.2,z=Z)
out5$p.value
out6$p.value
out7$p.value
#### Also can rerun with normal data to make sure continuous data work
set.seed(2)
Nc<-100
Nt<-200
Yt<-rnorm(Nt,mean=.3)
Yc<-rnorm(Nc)
Y<-c(Yt,Yc)
Z<-c(rep(2,Nt),rep(1,Nc))


out5b<-nicqTest(x=Y,g=nimDiffOR,delta0=.1,q=.2,z=Z)
out6b<-nicqTest(x=Y,g=nimDiff,delta0=.1,q=.2,z=Z)
out7b<-nicqTest(x=Y,g=nimOR,delta0=.1,q=.2,z=Z)

out5b$p.value
out6b$p.value
out7b$p.value


