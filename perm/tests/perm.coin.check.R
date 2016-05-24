## check perm against coin (ver 1.0-8 or greater)
library(perm)
library(coin)
packageDescription("coin")$Version

set.seed(1)
nsim<-3
#nsim<-100
n<-8
outperm1<-outcoin1<-rep(NA,nsim)
g<-factor(c(rep("a",n/2),rep("b",n/2)))
y<-rnorm(n)
permTS(y~g,method="exact.network",alternative="less")
independence_test(y~g,distribution=exact(),alternative="less")

permTS(y~g,method="pclt",alternative="less")
independence_test(y~g,distribution=asymptotic(),alternative="less")


permTS(y~g,method="exact.mc",alternative="less")
independence_test(y~g,distribution=approximate(),alternative="less")


for (i in 1:nsim){
    #y<-rpois(n,10)
    y<-rnorm(n)
    outperm1[i]<-permTS(y~g,method="exact.network",alternative="less")$p.value
    outcoin1[i]<-pvalue(independence_test(y~g,distribution=exact(),alternative="less"))
}

all.equal(outperm1,outcoin1)

n<-100
g<-factor(c(rep("a",n/2),rep("b",n/2)))
for (i in 1:nsim){
    #y<-rpois(n,10)
    y<-rnorm(n)
    outperm1[i]<-permTS(y~g,method="pclt",alternative="less")$p.value
    outcoin1[i]<-pvalue(independence_test(y~g,distribution=asymptotic(),alternative="less"))
}
all.equal(outperm1,outcoin1)

outperm2<-outcoin2<-rep(NA,nsim)
n<-100
set.seed(10301)
g<-factor(c(rep("a",n/4),rep("b",n/4),rep("c",n/4),rep("d",n/4)),levels=c("c","d","b","a"))
for (i in 1:nsim){
    #y<-rpois(n,10)
    y<-rnorm(n)
    outperm2[i]<-permKS(y~g,method="pclt")$p.value
    outcoin2[i]<-pvalue(independence_test(y~g,teststat="quad",distribution=asymptotic()))
}
all.equal(outperm2,outcoin2)


