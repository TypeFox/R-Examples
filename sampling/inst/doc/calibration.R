### R code from vignette source 'calibration.Snw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: calibration.Snw:20-23
###################################################
library(sampling)
ps.options(pointsize=12)
options(width=60)


###################################################
### code chunk number 2: calib1
###################################################
data = rbind(matrix(rep("A", 150), 150, 1, byrow = TRUE), 
matrix(rep("B", 100), 100, 1, byrow = TRUE))
data = cbind.data.frame(data, c(rep(1, 60), rep(2,50), rep(3, 60), rep(1, 40), rep(2, 40)),
 1000 * runif(250))
sex = runif(nrow(data))
for (i in 1:length(sex)) if (sex[i] < 0.3) sex[i] = 1 else sex[i] = 2
data = cbind.data.frame(data, sex)
names(data) = c("state", "region", "income", "sex")
summary(data)


###################################################
### code chunk number 3: calib2
###################################################
table(data$state)


###################################################
### code chunk number 4: calib3
###################################################
s=strata(data,c("state"),size=c(25,20), method="srswor")


###################################################
### code chunk number 5: calib31
###################################################
s=getdata(data,s)


###################################################
### code chunk number 6: calib32
###################################################
status=runif(nrow(s))
for(i in 1:length(status))
 if(status[i]<0.3) status[i]=0 else status[i]=1
s=cbind.data.frame(s,status) 


###################################################
### code chunk number 7: calib4
###################################################
s=rhg_strata(s,selection="region")


###################################################
### code chunk number 8: calib5
###################################################
sr=s[s$status==1,]


###################################################
### code chunk number 9: calib6
###################################################
X=cbind(disjunctive(data$sex),disjunctive(data$region))


###################################################
### code chunk number 10: calib7
###################################################
total=c(t(rep(1,nrow(data)))%*%X)


###################################################
### code chunk number 11: calib8
###################################################
Xs = X[sr$ID_unit,]
d = 1/(sr$Prob * sr$prob_resp)
summary(d)


###################################################
### code chunk number 12: calib9
###################################################
g = calib(Xs, d, total, method = "linear")
summary(g)


###################################################
### code chunk number 13: calib10
###################################################
w=d*g
summary(w)


###################################################
### code chunk number 14: calib11
###################################################
checkcalibration(Xs, d, total, g)


###################################################
### code chunk number 15: calib12
###################################################
cat("stratum 1\n")
data1=data[data$state=='A',]
X1=X[data$state=='A',]
total1=c(t(rep(1, nrow(data1))) %*% X1)
sr1=sr[sr$Stratum==1,]
Xs1=X[sr1$ID_unit,]
d1 = 1/(sr1$Prob * sr1$prob_resp)
g1=calib(Xs1, d1, total1, method = "linear")
checkcalibration(Xs1, d1, total1, g1)
cat("stratum 2\n")
data2=data[data$state=='B',]
X2=X[data$state=='B',]
total2=c(t(rep(1, nrow(data2))) %*% X2)
sr2=sr[sr$Stratum==2,]
Xs2=X[sr2$ID_unit,]
d2 = 1/(sr2$Prob * sr2$prob_resp)
g2=calib(Xs2, d2, total2, method = "linear")
checkcalibration(Xs2, d2, total2, g2)


###################################################
### code chunk number 16: calibration.Snw:156-175 (eval = FALSE)
###################################################
## data = rbind(matrix(rep("A", 150), 150, 1, byrow = TRUE), 
## matrix(rep("B", 100), 100, 1, byrow = TRUE))
## data = cbind.data.frame(data, c(rep(1, 60), rep(2,50), rep(3, 60), rep(1, 40), rep(2, 40)),
##  1000 * runif(250))
## sex = runif(nrow(data))
## for (i in 1:length(sex)) if (sex[i] < 0.3) sex[i] = 1 else sex[i] = 2
## data = cbind.data.frame(data, sex)
## names(data) = c("state", "region", "income", "sex")
## summary(data)
## table(data$state)
## s=strata(data,c("state"),size=c(25,20), method="srswor")
## s=getdata(data,s)
## status=runif(nrow(s))
## for(i in 1:length(status))
##  if(status[i]<0.3) status[i]=0 else status[i]=1
## s=cbind.data.frame(s,status) 
## s=rhg_strata(s,selection="region")
## sr=s[s$status==1,]
## X=cbind(disjunctive(data$sex),disjunctive(data$region))
## total=c(t(rep(1,nrow(data)))%*%X)
## Xs = X[sr$ID_unit,]
## d = 1/(sr$Prob * sr$prob_resp)
## summary(d)
## g = calib(Xs, d, total, method = "linear")
## summary(g)
## w=d*g
## summary(w)
## checkcalibration(Xs, d, total, g)
## cat("stratum 1\n")
## data1=data[data$state=='A',]
## X1=X[data$state=='A',]
## total1=c(t(rep(1, nrow(data1))) %*% X1)
## sr1=sr[sr$Stratum==1,]
## Xs1=X[sr1$ID_unit,]
## d1 = 1/(sr1$Prob * sr1$prob_resp)
## g1=calib(Xs1, d1, total1, method = "linear")
## checkcalibration(Xs1, d1, total1, g1)
## cat("stratum 2\n")
## data2=data[data$state=='B',]
## X2=X[data$state=='B',]
## total2=c(t(rep(1, nrow(data2))) %*% X2)
## sr2=sr[sr$Stratum==2,]
## Xs2=X[sr2$ID_unit,]
## d2 = 1/(sr2$Prob * sr2$prob_resp)
## g2=calib(Xs2, d2, total2, method = "linear")
## checkcalibration(Xs2, d2, total2, g2)
## 
## 
## 
## sampling.newpage()
## 


###################################################
### code chunk number 17: ex1
###################################################
X=cbind(c(rep(1,50),rep(0,50)),c(rep(0,50),rep(1,50)),1:100)
# vector of population totals
total=apply(X,2,"sum")
Z=150:249 
# the variable of interest
Y=5*Z*(rnorm(100,0,sqrt(1/3))+apply(X,1,"sum"))
# inclusion probabilities
pik=inclusionprobabilities(Z,20)
# joint inclusion probabilities 
pikl=UPtillepi2(pik)
# number of simulations; let nsim=10000 for an accurate result
nsim=10
c1=c2=c3=c4=c5=c6=numeric(nsim)
for(i in 1:nsim)
{
# draws a sample
s=UPtille(pik)
# computes the inclusion prob. for the sample
piks=pik[s==1]
# the sample matrix of auxiliary information
Xs=X[s==1,]
# computes the g-weights
g=calib(Xs,d=1/piks,total,method="linear")
# computes the variable of interest in the sample
Ys=Y[s==1]
# computes the joint inclusion prob. for the sample
pikls=pikl[s==1,s==1]
# computes the calibration estimator and its variance estimation
cc=calibev(Ys,Xs,total,pikls,d=1/piks,g,with=FALSE,EPS=1e-6)
c1[i]=cc$calest
c2[i]=cc$evar
# computes the variance estimator of the calibration estimator (second method)
c3[i]=varest(Ys,Xs,pik=piks,w=g/piks)
# computes the variance estimator of the HT estimator using varest()
c4[i]=varest(Ys,pik=piks) 
# computes the variance estimator of the HT estimator using varHT()
c5[i]=varHT(Ys,pikls,2) 
# computes the Horvitz-Thompson estimator
c6[i]=HTestimator(Ys,piks)
}
cat("the population total:",sum(Y),"\n")
cat("the calibration estimator under simulations:", mean(c1),"\n")
N=length(Y)
delta=matrix(0,N,N)
for(k in 1:(N-1))
  for(l in (k+1):N)
   delta[k,l]=delta[l,k]=pikl[k,l]-pik[k]*pik[l]
diag(delta)=pik*(1-pik)
var_HT=0
var_asym=0
e=lm(Y~X)$resid
for(k in 1:N)
  for(l in 1:N) {var_HT=var_HT+Y[k]*Y[l]*delta[k,l]/(pik[k]*pik[l])
                 var_asym=var_asym+e[k]*e[l]*delta[k,l]/(pik[k]*pik[l])} 
cat("the approximate variance of the calibration estimator:",var_asym,"\n")
cat("first variance estimator of the calibration est. using calibev function:\n")
cat("MSE of the first variance estimator:", var(c2)+(mean(c2)-var_asym)^2,"\n")
cat("second variance estimator of the calibration est. using varest function:\n")
cat("MSE of the second variance estimator:", var(c3)+(mean(c3)-var_asym)^2,"\n")
cat("the Horvitz-Thompson estimator under simulations:", mean(c6),"\n")
cat("the variance of the HT estimator:", var_HT, "\n")
cat("the variance estimator of the HT estimator under simulations:", mean(c4),"\n")
cat("MSE of the variance estimator 1 of HT estimator:", var(c4)+(mean(c4)-var_HT)^2,"\n")
cat("MSE of the variance estimator 2 of HT estimator:", var(c5)+(mean(c5)-var_HT)^2,"\n")


###################################################
### code chunk number 18: calibration.Snw:262-266 (eval = FALSE)
###################################################
## X=cbind(c(rep(1,50),rep(0,50)),c(rep(0,50),rep(1,50)),1:100)
## # vector of population totals
## total=apply(X,2,"sum")
## Z=150:249 
## # the variable of interest
## Y=5*Z*(rnorm(100,0,sqrt(1/3))+apply(X,1,"sum"))
## # inclusion probabilities
## pik=inclusionprobabilities(Z,20)
## # joint inclusion probabilities 
## pikl=UPtillepi2(pik)
## # number of simulations; let nsim=10000 for an accurate result
## nsim=10
## c1=c2=c3=c4=c5=c6=numeric(nsim)
## for(i in 1:nsim)
## {
## # draws a sample
## s=UPtille(pik)
## # computes the inclusion prob. for the sample
## piks=pik[s==1]
## # the sample matrix of auxiliary information
## Xs=X[s==1,]
## # computes the g-weights
## g=calib(Xs,d=1/piks,total,method="linear")
## # computes the variable of interest in the sample
## Ys=Y[s==1]
## # computes the joint inclusion prob. for the sample
## pikls=pikl[s==1,s==1]
## # computes the calibration estimator and its variance estimation
## cc=calibev(Ys,Xs,total,pikls,d=1/piks,g,with=FALSE,EPS=1e-6)
## c1[i]=cc$calest
## c2[i]=cc$evar
## # computes the variance estimator of the calibration estimator (second method)
## c3[i]=varest(Ys,Xs,pik=piks,w=g/piks)
## # computes the variance estimator of the HT estimator using varest()
## c4[i]=varest(Ys,pik=piks) 
## # computes the variance estimator of the HT estimator using varHT()
## c5[i]=varHT(Ys,pikls,2) 
## # computes the Horvitz-Thompson estimator
## c6[i]=HTestimator(Ys,piks)
## }
## cat("the population total:",sum(Y),"\n")
## cat("the calibration estimator under simulations:", mean(c1),"\n")
## N=length(Y)
## delta=matrix(0,N,N)
## for(k in 1:(N-1))
##   for(l in (k+1):N)
##    delta[k,l]=delta[l,k]=pikl[k,l]-pik[k]*pik[l]
## diag(delta)=pik*(1-pik)
## var_HT=0
## var_asym=0
## e=lm(Y~X)$resid
## for(k in 1:N)
##   for(l in 1:N) {var_HT=var_HT+Y[k]*Y[l]*delta[k,l]/(pik[k]*pik[l])
##                  var_asym=var_asym+e[k]*e[l]*delta[k,l]/(pik[k]*pik[l])} 
## cat("the approximate variance of the calibration estimator:",var_asym,"\n")
## cat("first variance estimator of the calibration est. using calibev function:\n")
## cat("MSE of the first variance estimator:", var(c2)+(mean(c2)-var_asym)^2,"\n")
## cat("second variance estimator of the calibration est. using varest function:\n")
## cat("MSE of the second variance estimator:", var(c3)+(mean(c3)-var_asym)^2,"\n")
## cat("the Horvitz-Thompson estimator under simulations:", mean(c6),"\n")
## cat("the variance of the HT estimator:", var_HT, "\n")
## cat("the variance estimator of the HT estimator under simulations:", mean(c4),"\n")
## cat("MSE of the variance estimator 1 of HT estimator:", var(c4)+(mean(c4)-var_HT)^2,"\n")
## cat("MSE of the variance estimator 2 of HT estimator:", var(c5)+(mean(c5)-var_HT)^2,"\n")
## 
## sampling.newpage()
## 


###################################################
### code chunk number 19: gen1
###################################################
N=400
n=100
X=rgamma(N,3,4)
total=sum(X)
Z=2*X+runif(N)
Y=3*X+rnorm(N)
print(cor(X,Y))
print(cor(X,Z))
L=1 
U=5
C=1.5  
A=(U-L)/((C-L)*(U-C))
p=((U-C)+(C-L)*exp(A*Y*0.3))/(L*(U-C)+U*(C-L)*exp(A*Y*0.3))
summary(p)
bounds=c(L,U)
s=srswor(n,N)
r=numeric(n)
for(j in 1:n) if(runif(1)<p[s==1][j]) r[j]=1  
print("Size of r is:")
nr=sum(r)
print(nr)
Xr=X[s==1][r==1]
Yr=Y[s==1][r==1]
Zr=Z[s==1][r==1]
pikr=rep(n/N,times=nr)
d=1/(pikr)
g1=gencalib(Xr,Zr,d,total,method="logit",bounds=bounds,C=C)
g2=gencalib(Xr,Yr,d,total,method="logit",bounds=bounds,C=C)
g3=gencalib(Xr,Xr,d,total,method="logit",bounds=bounds,C=C)
if(is.null(g1))
print("g1 is null") else
if(checkcalibration(Xr,d,total,g1)$result)
{print("the gen.calibration estimator using Zs as instrumental variable")
print(sum(Yr*g1*d))
}
if(is.null(g2))
print("g2 is null") else
if(checkcalibration(Xr,d,total,g2)$result)
{
print("the gen.calibration estimator using Ys as instrumental variable")
print(sum(Yr*g2*d))
}
if(is.null(g3))
print("g3 is null") else
if(checkcalibration(Xr,d,total,g3)$result)
{
print("the calibration estimator")
print(sum(Yr*g3*d))
}
cat("The population total is:", sum(Y),"\n")


###################################################
### code chunk number 20: calibration.Snw:390-394 (eval = FALSE)
###################################################
## N=400
## n=100
## X=rgamma(N,3,4)
## total=sum(X)
## Z=2*X+runif(N)
## Y=3*X+rnorm(N)
## print(cor(X,Y))
## print(cor(X,Z))
## L=1 
## U=5
## C=1.5  
## A=(U-L)/((C-L)*(U-C))
## p=((U-C)+(C-L)*exp(A*Y*0.3))/(L*(U-C)+U*(C-L)*exp(A*Y*0.3))
## summary(p)
## bounds=c(L,U)
## s=srswor(n,N)
## r=numeric(n)
## for(j in 1:n) if(runif(1)<p[s==1][j]) r[j]=1  
## print("Size of r is:")
## nr=sum(r)
## print(nr)
## Xr=X[s==1][r==1]
## Yr=Y[s==1][r==1]
## Zr=Z[s==1][r==1]
## pikr=rep(n/N,times=nr)
## d=1/(pikr)
## g1=gencalib(Xr,Zr,d,total,method="logit",bounds=bounds,C=C)
## g2=gencalib(Xr,Yr,d,total,method="logit",bounds=bounds,C=C)
## g3=gencalib(Xr,Xr,d,total,method="logit",bounds=bounds,C=C)
## if(is.null(g1))
## print("g1 is null") else
## if(checkcalibration(Xr,d,total,g1)$result)
## {print("the gen.calibration estimator using Zs as instrumental variable")
## print(sum(Yr*g1*d))
## }
## if(is.null(g2))
## print("g2 is null") else
## if(checkcalibration(Xr,d,total,g2)$result)
## {
## print("the gen.calibration estimator using Ys as instrumental variable")
## print(sum(Yr*g2*d))
## }
## if(is.null(g3))
## print("g3 is null") else
## if(checkcalibration(Xr,d,total,g3)$result)
## {
## print("the calibration estimator")
## print(sum(Yr*g3*d))
## }
## cat("The population total is:", sum(Y),"\n")
##                                  
## sampling.newpage()               
##                                  


