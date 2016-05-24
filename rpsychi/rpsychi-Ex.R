pkgname <- "rpsychi"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('rpsychi')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("FNONCT")
### * FNONCT

flush(stderr()); flush(stdout())

### Name: FNONCT
### Title: Compute the noncentrality parameter of an F distribution
### Aliases: FNONCT
### Keywords: internal

### ** Examples

##Kline (2004) Table 4.7
FNONCT(9.61, 1, 58, .975)
FNONCT(9.61, 1, 58, .025)



cleanEx()
nameEx("dep.oneway")
### * dep.oneway

flush(stderr()); flush(stdout())

### Name: dep.oneway
### Title: A one-way design with dependent samples using individual data:
###   Reporting effect size
### Aliases: dep.oneway
### Keywords: design htest

### ** Examples

##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=3)
                  )
dep.oneway(formula = y~x, data=dat, block="subj")


##contrast 1: a - c, contrast 2 : 1/2(a + c) - b
my.cont <- matrix(c(1,0,-1,1/2,-1,1/2), ncol=3, nrow=2, byrow=TRUE)
dep.oneway(formula = y~x, data=dat, block="subj", contr=my.cont)



cleanEx()
nameEx("dep.oneway.second")
### * dep.oneway.second

flush(stderr()); flush(stdout())

### Name: dep.oneway.second
### Title: A one-way design with dependent samples using published work:
###   Reporting effect size
### Aliases: dep.oneway.second
### Keywords: design htest

### ** Examples
 
##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=3)
                  )
dep.oneway(formula = y~x, data=dat, block="subj")

datwide <- reshape(dat, direction="wide", idvar="subj", timevar="x")
tmp <- datwide[,-1]
dep.oneway.second(m = apply(tmp, 2, mean), apply(tmp, 2, sd), n = nrow(tmp), corr=cor(tmp))



##Kline (2004) Table 6.15
my.cont <- matrix(c(-5,-3,-1,1,3,5,
                   5,-1,-4,-4,-1,5), ncol=6, nrow=2, byrow=TRUE)
dep.oneway.second(m = c(11.77,21.39,27.5,31.02,32.58,34.2), 
                  sd = c(7.6,8.44,8.95,9.21,9.49,9.62), 
                  n = 137, 
                  corr=lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
                  .84,.8,.91,.88,.93)),
                  contr=my.cont)



cleanEx()
nameEx("dep.t.test")
### * dep.t.test

flush(stderr()); flush(stdout())

### Name: dep.t.test
### Title: A t-test with dependent samples using individual data: Reporting
###   effect size
### Aliases: dep.t.test
### Keywords: design htest

### ** Examples
 
##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=2)
                  )
dep.t.test(y~x, block="subj", data=dat)




cleanEx()
nameEx("dep.t.test.second")
### * dep.t.test.second

flush(stderr()); flush(stdout())

### Name: dep.t.test.second
### Title: A t-test with dependent samples using published work: Reporting
###   effect size
### Aliases: dep.t.test.second
### Keywords: design htest

### ** Examples
 
##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5),
                  subj = rep(paste("s", 1:5, sep=""), times=2)
                  )
datwide <- reshape(dat, direction="wide", idvar="subj", timevar="x")

dep.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = nlevels(dat$subj),
                  corr = cor(datwide[,2:3])[1,2]
                  )

dep.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = 30,
                  corr = cor(datwide[,2:3])[1,2]
                  )



cleanEx()
nameEx("formatted")
### * formatted

flush(stderr()); flush(stdout())

### Name: formatted
### Title: Convert a numeric vector into a character vector with the
###   specified number of decimal place
### Aliases: formatted
### Keywords: univar

### ** Examples

data(infert)
x <- svar(infert$age)    #sample variance
formatted(x)
formatted(x, digits=4)



cleanEx()
nameEx("groupSummary")
### * groupSummary

flush(stderr()); flush(stdout())

### Name: groupSummary
### Title: Compute summary statistics by group
### Aliases: groupSummary
### Keywords: univar

### ** Examples

data(infert)
infert$case <- factor(infert$case, labels=c("control", "case"))
infert$induced <- factor(infert$induced, labels=c("0","1","2 or more"))
infert$spontaneous <- factor(infert$spontaneous, labels=c("0","1","2 or more"))

#continuous and categorical variables
groupSummary(infert, group="case")

#continuous variables only
groupSummary(infert[, c(2,3,7,8, 5)], group="case")

#categorical variables only
groupSummary(infert[, c(1,4, 6, 5)],  group="case")    

#total sample
groupSummary(infert[, c(1,4, 6, 5)])




cleanEx()
nameEx("ind.oneway")
### * ind.oneway

flush(stderr()); flush(stdout())

### Name: ind.oneway
### Title: A one-way design with independent samples using individual data:
###   Reporting effect size
### Aliases: ind.oneway
### Keywords: design htest

### ** Examples

##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5)
                  )                 
ind.oneway(formula = y~x, data=dat, sig.level=.05, digits=3)


##contrast 1: a - c, contrast 2: 1/2(a + c) - b
my.cont <- matrix(c(1,0,-1,1/2,-1,1/2), ncol=3, nrow=2, byrow=TRUE)
ind.oneway(formula = y~x, data=dat, contr=my.cont, sig.level=.05, digits=3)




cleanEx()
nameEx("ind.oneway.second")
### * ind.oneway.second

flush(stderr()); flush(stdout())

### Name: ind.oneway.second
### Title: A one-way design with independent samples using published work:
###   Reporting effect size
### Aliases: ind.oneway.second
### Keywords: design htest

### ** Examples

##Kline (2004) Table 6.3
dat <- data.frame(y = c(9,12,13,15,16,
                       8,12,11,10,14,
                       10,11,13,11,15),
                  x =  rep(factor(c("a","b","c")), each=5)
                  )                 

##contrast 1: a - c, contrast 2: 1/2(a + c) - b
my.cont <- matrix(c(1,0,-1,1/2,-1,1/2), ncol=3, nrow=2, byrow=TRUE)


ind.oneway.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n= tapply(dat$y, dat$x, length)) 

ind.oneway.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n= tapply(dat$y, dat$x, length),
                  contr = my.cont)   




cleanEx()
nameEx("ind.prop")
### * ind.prop

flush(stderr()); flush(stdout())

### Name: ind.prop
### Title: A Z test for the equality of two proportions using individual
###   data: Reporting effect size
### Aliases: ind.prop
### Keywords: design htest

### ** Examples

##Kline (2004) Chapter 5
x1 <- c("relapsed", "not relapsed")
y1 <- c("control", "treatment")

dat <- data.frame(y =         
factor(c(rep(x1, c(60, 40)), rep(x1, c(40, 60))), levels=x1),
x = factor(rep(y1, each=100), levels=y1)
)
tab <- xtabs(~x+y, data=dat)
tab
ind.prop(y~x, data=dat, lev.count=2, ref.ind=1)    #Odds for not relapse is higher in treatment than control condition.
ind.prop(y~x, data=dat, lev.count=1, ref.ind=1)    #Odds for relapse is lower in treatment than control condition.
ind.prop(y~x, data=dat, lev.count=2, ref.ind=2)    #Odds for not relapse is lower in control than treatment condition.
ind.prop(y~x, data=dat, lev.count=1, ref.ind=2)    #Odds for relapse is higher in control than treatment condition.



cleanEx()
nameEx("ind.prop.second")
### * ind.prop.second

flush(stderr()); flush(stdout())

### Name: ind.prop.second
### Title: A Z test for the equality of two proportions using published
###   workc
### Aliases: ind.prop.second
### Keywords: design htest

### ** Examples

##Kline (2004) Chapter 5
x1 <- c("relapsed", "not relapsed")
y1 <- c("control", "treatment")

dat <- data.frame(y =         
factor(c(rep(x1, c(60, 40)), rep(x1, c(40, 60))), levels=x1),
x = factor(rep(y1, each=100), levels=y1)
)
tab <- xtabs(~x+y, data=dat)
tab
ind.prop.second(x=tab[,1], n = rowSums(tab))             #Risk for relapse is lower in treatment than control condition.
ind.prop.second(x=tab[,1], n = rowSums(tab), ref.ind=2)  #Risk for relapse is higher in control than treatment condition.



cleanEx()
nameEx("ind.t.test")
### * ind.t.test

flush(stderr()); flush(stdout())

### Name: ind.t.test
### Title: A t-test with independent samples using individual data:
###   Reporting effect size
### Aliases: ind.t.test
### Keywords: design htest

### ** Examples

##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5)
                  )
ind.t.test(y~x, data=dat, correct=FALSE)



cleanEx()
nameEx("ind.t.test.second")
### * ind.t.test.second

flush(stderr()); flush(stdout())

### Name: ind.t.test.second
### Title: A t-test with independent samples using published work:
###   Reporting effect size
### Aliases: ind.t.test.second
### Keywords: design htest

### ** Examples

##Kline (2004) Table 4.4
dat <- data.frame(y = c(9,12,13,15,16,8,12,11,10,14),
                  x =  rep(factor(c("a","b")), each=5)
                  )
ind.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = tapply(dat$y, dat$x, length), correct=FALSE
                  )
ind.t.test.second(m = tapply(dat$y, dat$x, mean),
                  sd = tapply(dat$y, dat$x, sd),
                  n = tapply(dat$y, dat$x, length), correct=TRUE
                  )     #approximate unbiased estimator of delta



cleanEx()
nameEx("ind.twoway")
### * ind.twoway

flush(stderr()); flush(stdout())

### Name: ind.twoway
### Title: A two-way design with independent samples using individual data
### Aliases: ind.twoway
### Keywords: design htest

### ** Examples

##Kline (2004) Table 7.5
dat <- data.frame(
           y = c(2,3,4,1,3,1,3,4,5,5,6,6,6,7),
           A = factor(c(rep("A1",5), rep("A2", 9))),
           B = factor(c(rep("B1",3), rep("B2",2), rep("B1",2), rep("B2",7)))
           )

ind.twoway(y~A*B, data=dat)



cleanEx()
nameEx("ind.twoway.second")
### * ind.twoway.second

flush(stderr()); flush(stdout())

### Name: ind.twoway.second
### Title: A two-way design with independent samples using published work
### Aliases: ind.twoway.second
### Keywords: design htest

### ** Examples

##Cohen (2000) Table 1
m.mat  <- matrix(c(37.13, 39.31, 39.22, 32.71), ncol=2) #2 * 2
sd.mat <- matrix(c(13.82, 9.42, 9.43, 9.62), ncol=2)
n.mat <- matrix(c(9, 13, 8, 14), ncol=2)

ind.twoway.second(m = m.mat, sd = sd.mat, n = n.mat)


##Tabachnick and Fidell (2007) 
#5.7 Complete example of two-way randomized-groups ANOVA (p.221-236)
m.mat <- matrix(c(837.9, 573.6, 354.9, 699.0, 112.0, 
      852.2, 781.6, 683.3, 1193.9, 130.0), ncol=2)    #5 * 2
sd.mat <- matrix(c(189.87449, 61.31195, 147.93351, 128.51891, 43.36922, 
      227.17042, 104.81221, 116.25934, 198.36692, 37.64158), ncol=2) #5 * 2
n.mat <- matrix(rep(10, 10), ncol=2)

ind.twoway.second(m = m.mat, sd = sd.mat, n = n.mat)


##Kline (2004) Table 7.5
dat <- data.frame(
           y = c(2,3,4,1,3,1,3,4,5,5,6,6,6,7),
           A = factor(c(rep("A1",5), rep("A2", 9))),
           B = factor(c(rep("B1",3), rep("B2",2), rep("B1",2), rep("B2",7)))
           )
ind.twoway.second(m = tapply(dat$y, list(dat$A,dat$B), mean), 
                  sd = tapply(dat$y, list(dat$A,dat$B), sd), 
                  n = tapply(dat$y, list(dat$A,dat$B), length)
                    )



cleanEx()
nameEx("lower2R")
### * lower2R

flush(stderr()); flush(stdout())

### Name: lower2R
### Title: Convert a vector containing correlations into a correlation
###   matrix
### Aliases: lower2R
### Keywords: array

### ** Examples

lower2R(c(1:15))

##Kline (2004) Table 6.15
lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
        .84,.8,.91,.88,.93))

lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
        .84,.8,.91,.88,.93), 
        varname=paste("trial", 1:6, sep=""))



cleanEx()
nameEx("md.pattern2")
### * md.pattern2

flush(stderr()); flush(stdout())

### Name: md.pattern2
### Title: Display missing-data patterns
### Aliases: md.pattern2
### Keywords: array

### ** Examples

##Iwasaki (2002)
dat <- data.frame(matrix(c(
71  ,  68,  72,  72,  90,  72,  77,  76,  84,  77,
1850,2000,2100,1700,  NA,2200,2150,  NA,  NA,  NA,
136 , 139, 147, 142,  NA, 150, 156,  NA, 152,  NA,
34  , 45 ,  50,  38,  NA,  41,  43,  52,  57,  48
), ncol=4))
md.pattern2(dat)


#sample sizes in the specific pattern
#^
#^               numbers of missing data in each pattern
#|                ^
#|                |
#    X2 X3 X4 X1 NA
#6    1  1  1  1  0
#2    0  0  1  1  2
#1    0  0  0  1  3
#1    0  1  1  1  1
#Sum  4  3  1  0  8 --> numbers of missing data in each variable




cleanEx()
nameEx("multreg")
### * multreg

flush(stderr()); flush(stdout())

### Name: multreg
### Title: A multiple regression analysis using individual data
### Aliases: multreg
### Keywords: design htest

### ** Examples

##Cohen (2003) Table 3.5.1
dat <- data.frame(
salary = c(51876, 54511, 53425, 61863, 52926, 47034, 66432, 61100, 41934, 
  47454, 49832, 47047, 39115, 59677, 61458, 54528, 60327, 56600, 
  52542, 50455, 51647, 62895, 53740, 75822, 56596, 55682, 62091, 
  42162, 52646, 74199, 50729, 70011, 37939, 39652, 68987, 55579, 
  54671, 57704, 44045, 51122, 47082, 60009, 58632, 38340, 71219, 
  53712, 54782, 83503, 47212, 52840, 53650, 50931, 66784, 49751, 
  74343, 57710, 52676, 41195, 45662, 47606, 44301, 58582),
pubs  = c(18, 3, 2, 17, 11, 6, 38, 48, 9, 22, 30, 21, 
  10, 27, 37, 8, 13, 6, 12, 29, 29, 7, 6, 69, 11, 9, 
  20, 41, 3, 27, 14, 23, 1, 7, 19, 11, 31, 9, 12, 32, 
  26, 12, 9, 6, 39, 16, 12, 50, 18, 16, 5, 20, 50, 
  6, 19, 11, 13, 3, 8, 11, 25, 4),
cits = c(50, 26, 50, 34, 41, 37, 48, 56, 19, 29, 
    28, 31, 25, 40, 61, 32, 36, 69, 47, 29, 35, 
    35, 18, 90, 60, 30, 27, 35, 14, 56, 50, 25, 
    35, 1, 69, 69, 27, 50, 32, 33, 45, 54, 47, 29, 
    69, 47, 43, 55, 33, 28, 42, 24, 31, 27, 
    83, 49, 14, 36, 34, 70, 27, 28)   
)
multreg(salary~ pubs + cits, data=dat)



cleanEx()
nameEx("multreg.second")
### * multreg.second

flush(stderr()); flush(stdout())

### Name: multreg.second
### Title: A multiple regression analysis using published work
### Aliases: multreg.second
### Keywords: design htest

### ** Examples

##Cohen (2003) Table 3.5.1
dat <- data.frame(
salary = c(51876, 54511, 53425, 61863, 52926, 47034, 66432, 61100, 41934, 
  47454, 49832, 47047, 39115, 59677, 61458, 54528, 60327, 56600, 
  52542, 50455, 51647, 62895, 53740, 75822, 56596, 55682, 62091, 
  42162, 52646, 74199, 50729, 70011, 37939, 39652, 68987, 55579, 
  54671, 57704, 44045, 51122, 47082, 60009, 58632, 38340, 71219, 
  53712, 54782, 83503, 47212, 52840, 53650, 50931, 66784, 49751, 
  74343, 57710, 52676, 41195, 45662, 47606, 44301, 58582),
pubs  = c(18, 3, 2, 17, 11, 6, 38, 48, 9, 22, 30, 21, 
  10, 27, 37, 8, 13, 6, 12, 29, 29, 7, 6, 69, 11, 9, 
  20, 41, 3, 27, 14, 23, 1, 7, 19, 11, 31, 9, 12, 32, 
  26, 12, 9, 6, 39, 16, 12, 50, 18, 16, 5, 20, 50, 
  6, 19, 11, 13, 3, 8, 11, 25, 4),
cits = c(50, 26, 50, 34, 41, 37, 48, 56, 19, 29, 
    28, 31, 25, 40, 61, 32, 36, 69, 47, 29, 35, 
    35, 18, 90, 60, 30, 27, 35, 14, 56, 50, 25, 
    35, 1, 69, 69, 27, 50, 32, 33, 45, 54, 47, 29, 
    69, 47, 43, 55, 33, 28, 42, 24, 31, 27, 
    83, 49, 14, 36, 34, 70, 27, 28)   )

multreg.second(salary~ pubs + cits, corr=cor(dat), n= nrow(dat))
multreg.second(salary~ pubs + cits, corr=cor(dat), n= nrow(dat), 
        m = apply(dat, 2, mean), sd=apply(dat, 2, sd))



cleanEx()
nameEx("power.f")
### * power.f

flush(stderr()); flush(stdout())

### Name: power.f
### Title: Statistical power for F tests on means in the analysis of
###   variance
### Aliases: power.f
### Keywords: internal

### ** Examples

##Cohen (1988) ex.8.1
power.f(u = 3, n = 20, delta=.28, sig.level=.05)
power.f(u = 3, n = 20, delta=.28, sig.level=.10)


##Cohen (1988) ex 8.2
power.f(u = 2, n = 200, delta=.23, sig.level=.01)
power.f(u = 2, n = 200, delta=.33, sig.level=.01)



cleanEx()
nameEx("power.f2")
### * power.f2

flush(stderr()); flush(stdout())

### Name: power.f2
### Title: Statistical power for F tests on means in the two-way analysis
###   of variance
### Aliases: power.f2
### Keywords: internal

### ** Examples

##Cohen (1988) ex 8.6
power.f2(df1=1,df2=96,sig.level=0.05,delta=.10)
power.f2(df1=2,df2=96,sig.level=0.05,delta=.25)
power.f2(df1=3,df2=96,sig.level=0.05,delta=.40)

power.f2(df1=1,df2=120,sig.level=0.05,delta=.10)
power.f2(df1=2,df2=120,sig.level=0.05,delta=.25)
power.f2(df1=3,df2=120,sig.level=0.05,delta=.40)



cleanEx()
nameEx("power.multi")
### * power.multi

flush(stderr()); flush(stdout())

### Name: power.multi
### Title: Statistical power for F tests of variance proportions
### Aliases: power.multi
### Keywords: internal

### ** Examples

##Cohen (1988) ex 9.1
power.multi(n=95, n.ind=5, delta=.1111, sig.level=.05)




cleanEx()
nameEx("power.prop")
### * power.prop

flush(stderr()); flush(stdout())

### Name: power.prop
### Title: Statistical power for differences between proportions
### Aliases: power.prop
### Keywords: internal

### ** Examples

##Cohen (1988) ex.6.2
power.prop(h=.4, n=c(100,100), sig.level=0.05)



cleanEx()
nameEx("power.r")
### * power.r

flush(stderr()); flush(stdout())

### Name: power.r
### Title: Statistical power for the significance testing of a product
###   moment correlation
### Aliases: power.r
### Keywords: internal

### ** Examples

##Cohen (1988) ex.3.1
power.r(n=50, delta=0.3, sig.level=0.05)



cleanEx()
nameEx("power.t")
### * power.t

flush(stderr()); flush(stdout())

### Name: power.t
### Title: Statistical power for the t test for means
### Aliases: power.t
### Keywords: internal

### ** Examples

##Cohen (1988) ex.2.1
power.t(sig.level=.05, delta=0.5, n1=30, n2=30)



cleanEx()
nameEx("r2cov")
### * r2cov

flush(stderr()); flush(stdout())

### Name: r2cov
### Title: Convert correlation matrix into covariance matrix
### Aliases: r2cov
### Keywords: array

### ** Examples
 
##data(iris) 
x <- iris[,1:4] 
cov(x)
r2cov(apply(x, 2, sd), cor(x)) 


##Toyoda (1998) p.34 
r2cov(sd = sqrt(c(.862, 1.089, 0.606)), 
      R = lower2R(c(.505, -0.077, -.233)))



cleanEx()
nameEx("rpsychi-package")
### * rpsychi-package

flush(stderr()); flush(stdout())

### Name: rpsychi-package
### Title: Statistics for psychiatric research
### Aliases: rpsychi-package rpsychi
### Keywords: package

### ** Examples

##Kline (2004) Table 6.15
my.cont <- matrix(c(-5,-3,-1,1,3,5,
                   5,-1,-4,-4,-1,5), ncol=6, nrow=2, byrow=TRUE)
dep.oneway.second(m = c(11.77,21.39,27.5,31.02,32.58,34.2), 
                  sd = c(7.6,8.44,8.95,9.21,9.49,9.62), 
                  n = 137, 
                  corr=lower2R(c(.77,.59,.50,.48,.46,.81,.72,.69,.68,.89,
                  .84,.8,.91,.88,.93)),
                  contr=my.cont)



cleanEx()
nameEx("samplesize.d")
### * samplesize.d

flush(stderr()); flush(stdout())

### Name: samplesize.d
### Title: Sample size estimation in the t test for means
### Aliases: samplesize.d
### Keywords: design htest

### ** Examples

##Cohen (1988) ex.2.9
samplesize.d(delta=.20, power=.95, sig.level=.05)



cleanEx()
nameEx("samplesize.etasq")
### * samplesize.etasq

flush(stderr()); flush(stdout())

### Name: samplesize.etasq
### Title: Sample size estimation in the analysis of variance
### Aliases: samplesize.etasq
### Keywords: design htest

### ** Examples

##Cohen (1988) ex.8.10
f <- .25
samplesize.etasq(k=4, delta= f^2/(1+f^2), power=.80, sig.level=.05)



cleanEx()
nameEx("samplesize.h")
### * samplesize.h

flush(stderr()); flush(stdout())

### Name: samplesize.h
### Title: Sample size estimation in the differences between proportions
### Aliases: samplesize.h
### Keywords: design htest

### ** Examples

##Cohen (1988) ex.6.7
samplesize.h(delta=.20, power=.90, sig.level=.01)

##Cohen (1988) ex.6.8
samplesize.h(delta=.2828, power=.95, sig.level=.05)



cleanEx()
nameEx("samplesize.r")
### * samplesize.r

flush(stderr()); flush(stdout())

### Name: samplesize.r
### Title: Sample size estimation in the significance testing of a product
###   moment correlation
### Aliases: samplesize.r
### Keywords: design htest

### ** Examples

##Cohen (1988) ex.3.4
samplesize.r(delta=.3, power=.8, sig.level=.05)



cleanEx()
nameEx("samplesize.rsq")
### * samplesize.rsq

flush(stderr()); flush(stdout())

### Name: samplesize.rsq
### Title: Sample size estimation in the F tests of variance proportions
### Aliases: samplesize.rsq
### Keywords: design htest

### ** Examples

##Cohen (1988) ex.9.18
samplesize.rsq(delta=.16, n.ind=20, power = .90, sig.level=.05)



cleanEx()
nameEx("ssd")
### * ssd

flush(stderr()); flush(stdout())

### Name: ssd
### Title: Compute sample standard deviation
### Aliases: ssd
### Keywords: univar

### ** Examples

data(infert)
ssd(infert$age)     #sample standard deviation
sd(infert$age)      #unbiased standard deviation



cleanEx()
nameEx("ssd2sd")
### * ssd2sd

flush(stderr()); flush(stdout())

### Name: ssd2sd
### Title: Convert sample standard deviation into unbiased one
### Aliases: ssd2sd
### Keywords: univar

### ** Examples

data(infert)
ssd2sd(nrow(infert), ssd(infert$age))
sd(infert$age)



cleanEx()
nameEx("svar")
### * svar

flush(stderr()); flush(stdout())

### Name: svar
### Title: Compute sample variance
### Aliases: svar
### Keywords: univar

### ** Examples

data(infert)
svar(infert$age)    #sample variance
var(infert$age)     #unbiased variance



cleanEx()
nameEx("zero.r.test")
### * zero.r.test

flush(stderr()); flush(stdout())

### Name: zero.r.test
### Title: A significance testing of a product moment correlation using
###   individual data
### Aliases: zero.r.test
### Keywords: design htest

### ** Examples

dat <- data.frame(x = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1),
                  y = c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8))          
zero.r.test(y~x, data=dat)



cleanEx()
nameEx("zero.r.test.second")
### * zero.r.test.second

flush(stderr()); flush(stdout())

### Name: zero.r.test.second
### Title: A significance testing of a product moment correlation using
###   published work
### Aliases: zero.r.test.second
### Keywords: design htest

### ** Examples

zero.r.test.second(r = 0.571, n = 9)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
