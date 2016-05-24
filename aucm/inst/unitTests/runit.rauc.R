library("RUnit")
library("aucm")

test.rauc <- function() {

RNGkind("Mersenne-Twister", "Inversion")

tol.=1e-6
# rauc call result depends on the state of rng. here we simulate data before each call to rauc to guarrantee reproducibility

dat = sim.dat.1(n=200,seed=1)
checkEqualsNumeric(mean(as.matrix(dat)), 0.1411382, tol=tol.)

# BYlogreg
fit=rlogit(y~x1, dat)
checkEqualsNumeric(coef(fit), c(0.008824577, 1.742154014), tol=tol.)
fit=rlogit(y~-1+x1, dat) # test intercept in formula
checkEqualsNumeric(coef(fit), c(0.008824577, 1.742154014), tol=tol.)

# sauc
fit = sauc.phi(y~x1+x2, dat, constrain.method="L2", h.method="Lin", start.method="rlogit", ret.vcov=F) 
checkEqualsNumeric(coef(fit), c(0.8468809, -0.5317826), tol=tol.)

# raucl
fit.1 = rauc (y~x1+x2, dat, lambda=2, kernel="linear", maxit=2, minQuad.control = control.minQuad(method = 'default', ws="v"), verbose=TRUE)
checkEqualsNumeric(coef(fit.1), c(4.691588, -2.964932), tolerance=tol.)
checkEqualsNumeric(predict(fit.1, dat)[1:3], c(-0.1437135, -0.1943497,  0.0770122), tolerance=tol.)

# raucr
fit.1 = rauc (y~x1+x2, dat, lambda=2, kernel="rbf", maxit=2, minQuad.control = control.minQuad(method = 'default', ws="v"), verbose=TRUE, para=1)
checkEqualsNumeric(trainauc(fit.1), 0.7393958, tolerance=tol.)
checkEqualsNumeric(predict(fit.1, dat)[1:3], c(-0.04663125, -0.45071539,  0.57591081), tolerance=tol.)


## quadractic kernel
#
#dat=sim.NL1(n=200, seed=1)
#dat.test=sim.NL1(n=1e4, seed=1001)
#
#fit.l = rauc (y~x1+x2+x3+x4, dat, lambda=2, kernel="linear", maxit=1e3)
#fast.auc(predict(fit.l, dat.test), dat.test$y)
#
#fit.p = rauc (y~x1+x2+x3+x4, dat, lambda=2, kernel="p", para=2, maxit=1e3)
#fast.auc(predict(fit.p, dat.test), dat.test$y)
#
#require(svmw)
#fit.s = svml (y~x1+x2+x3+x4, dat, fitted=T, kernel="l", cost=1)
#score = predict(fit.s, dat.test)[[2]][,1]
#fast.auc(score, dat.test$y)
#
#fit.p = svml (y~x1+x2+x3+x4, dat, fitted=F, kernel="p", degree=2, cost=1)
#score = predict(fit.p, dat.test)[[2]][,1]
#fast.auc(score, dat.test$y)
#


#
#
#fit1 = rauc (y~x1+x2, dat, lambda=1, kernel="linear", verbose=TRUE)
#fit1$training.auc # 0.7210025
#
#fit1 = rauc (y~x1+x2, dat, lambda=2, kernel="linear", verbose=TRUE, s=2)
#fit1$training.auc #  0.7207018
#
#
#
#fit3 = rauc (y~x1+x2, dat, lambda=2, kernel="rbf", para=1, verbose=TRUE)
#fit3$training.auc # 0.7773434
#
#fit4 = svml (y~x1+x2, dat, kernel="r", fitted=FALSE, cost=1e4) 
#fast.auc(predict(fit4, dat)$posterior[,1], dat$y) # 0.7921805
#tune.svml(y~x1+x2, dat, kernel="r")
##        1        10       100      1000     10000     1e+05
##0.7027569 0.7254135 0.7517794 0.7653133 0.7921805 0.6674687
#
## glm derived score for comparision
#fit.glm=glm(y~x1+x2, dat, family="binomial")
#fast.auc(fit1$X %*% fit.glm$coef[-1], fit1$y) # 
#
## add outliers
#dat = sim.dat.1(n=200,seed=1, add.outliers=TRUE)
#
#fit3 = rauc (y~x1+x2, dat, lambda=2, kernel="rbf", para=1, verbose=TRUE)
#fit3$training.auc # 0.7066667
#
#fit4 = svml (y~x1+x2, dat, kernel="r", fitted=FALSE, cost=1e4) 
#fast.auc(predict(fit4, dat)$posterior[,1], dat$y) # 0.6910101
#tune.svml(y~x1+x2, dat, kernel="r")
##        1        10       100      1000     10000     1e+05 
##0.6485859 0.6705051 0.6722222 0.6767677 0.6910101 0.5007071
#
#
## performance bench mark
#system.time(rauc (y~x1+x2, sim.dat.1(n=200,seed=1), lambda=2, kernel="linear", verbose=TRUE))
## output
##iter 1, beta[1] 0.9932739, beta[-1]/beta[1] -0.6899373, satured 0.4258647, val 0.3035524
##iter 2, beta[1] 1.403117, beta[-1]/beta[1] -0.6761924, satured 0.5804511, val 0.2906623
##iter 3, beta[1] 1.77679, beta[-1]/beta[1] -0.6732515, satured 0.6638596, val 0.2857926
##iter 4, beta[1] 2.151545, beta[-1]/beta[1] -0.6722024, satured 0.7181955, val 0.2830593
##iter 5, beta[1] 2.436216, beta[-1]/beta[1] -0.6680533, satured 0.7544862, val 0.281683
##iter 6, beta[1] 2.733543, beta[-1]/beta[1] -0.6655527, satured 0.7799499, val 0.2806743
##iter 7, beta[1] 2.931791, beta[-1]/beta[1] -0.6623021, satured 0.7957895, val 0.2801885
##iter 8, beta[1] 3.054115, beta[-1]/beta[1] -0.6603076, satured 0.8039098, val 0.2799963
##iter 9, beta[1] 3.122357, beta[-1]/beta[1] -0.6592986, satured 0.8081203, val 0.2799076
##iter 10, beta[1] 3.155542, beta[-1]/beta[1] -0.6593205, satured 0.8107268, val 0.2798711
##iter 11, beta[1] 3.173086, beta[-1]/beta[1] -0.6595918, satured 0.8122306, val 0.2798521
##iter 12, beta[1] 3.18549, beta[-1]/beta[1] -0.6598627, satured 0.8128321, val 0.2798396
##iter 13, beta[1] 3.19247, beta[-1]/beta[1] -0.6600143, satured 0.8133333, val 0.2798327
##[1] "total time used: 1.66856683095296 mins"
##   user  system elapsed 
##  99.30    0.23  100.24 
#
#
#
## larger n
#dat = sim.dat.1(n=300,seed=1)
#fit1 = rauc (y~x1+x2, dat, lambda=2, kernel="linear", verbose=TRUE)
#
#
############################################################
## a nonlinear example
#
#dat=skin.orange (n=100,seed=1,noise=FALSE)
#dim(dat)
#
## nonlinear kernel fit
#fit1 = rauc (y~x1+x2+x3+x4, dat, lambda=2, kernel="rbf", para=1, verbose=TRUE)
## glm fit
#fit.glm=glm(y~x1+x2+x3+x4, dat, family="binomial")
## linear kernel fit
#fit2 = rauc (y~x1+x2+x3+x4, dat, lambda=2, kernel="linear", start.method = "rlogit", verbose=TRUE)
#
## training data prediction
#fast.auc(fit1$linear.combination, fit1$y)
#fast.auc(fit1$X %*% fit.glm$coef[-1], fit1$y)
#fast.auc(fit2$linear.combination, fit2$y)
#
## test data prediction
#newdata=skin.orange (n=1000,seed=2,noise=FALSE)
#fast.auc(predict(fit1, newdata), newdata$y)
#fast.auc(as.matrix(subset(newdata, select=c(x1,x2,x3,x4))) %*% fit.glm$coef[-1], newdata$y)
#fast.auc(predict(fit2, newdata), newdata$y)
#
#
#
########################################################### IMPROVEMENTS #####################################################
#
# 
### rank = 2 problem 
#dat = sim.dat.1(n=300,seed=1,add.outliers = TRUE,std.dev = 1.0);fm = y~x1+x2
#
### linear kernel and random working set selection - low rank (2) problem
### setting initial alpha (to be passed to minQuad at each iteration in dca-loop) to estimate from previous dca() iteration 
### size of working set is automatically set
#set.seed(100) 
#fit.lin = rauc (fm, dat,lambda=.1,kernel="linear",
#verbose=TRUE,maxit = 100,tol = 1e-5,
#init.alpha.from.previous = TRUE,mem.efficient = TRUE,
#minQuad.control = control.minQuad(
#                            verbose = 1,maxit = 1e6,tol = 1e-4,
#                            method = "tron",                            
#                            working.set= "rv2wg")
#)
#
### 'rbf' kernel and random working set selection
### low rank mapped to possibly infinite rank problem try larger working set 'q' set.seed(100) 
### size of working set is set to q = 100
#fit.rbf = rauc (fm, dat,lambda=.1,kernel="rbf",para = 1, verbose=TRUE,maxit = 100,tol = 1e-5,
#init.alpha.from.previous = TRUE,mem.efficient = TRUE,
#minQuad.control = control.minQuad(
#                            verbose = 1,maxit = 1e6,tol = 1e-4,
#                            q = 100,
#                            method = "tron",                            
#                            working.set= "rv2wg")
#)


}
