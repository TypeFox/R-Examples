### R code from vignette source 'HT_Hajek_estimators.Snw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: HT_Hajek_estimators.Snw:22-26
###################################################
library(sampling)
ps.options(pointsize=12)
options(width=60)



###################################################
### code chunk number 2: up1
###################################################
data(belgianmunicipalities)
attach(belgianmunicipalities)
# sample size
n=20
pik=inclusionprobabilities(Tot04,n)
N=length(pik)




###################################################
### code chunk number 3: up2
###################################################
sim=10
ss=ss1=array(0,c(sim,4))



###################################################
### code chunk number 4: up3
###################################################
cat("Case 1\n")
y1=rep(3,N)
cat("Case 2\n")
y2=TaxableIncome
cat("Case 3\n")
x=1:N
pik3=inclusionprobabilities(x,n)
y3=1/pik3
cat("Case 4\n")
epsilon=rnorm(N,0,sqrt(1/3))
pik4=pik3
y4=5*(x+epsilon)



###################################################
### code chunk number 5: up4
###################################################
ht=numeric(4)
hajek=numeric(4)
for(i in 1:sim)
{
cat("Simulation ",i,"\n")
cat("Case 1\n")
s=UPtille(pik)
ht[1]=HTestimator(y1[s==1],pik[s==1])
hajek[1]=Hajekestimator(y1[s==1],pik[s==1],N,type="total")
cat("Case 2\n")
s1=UPpoisson(pik)
ht[2]=HTestimator(y2[s1==1],pik[s1==1])
hajek[2]=Hajekestimator(y2[s1==1],pik[s1==1],N,type="total")
cat("Case 3\n")
ht[3]=HTestimator(y3[s==1],pik3[s==1])
hajek[3]=Hajekestimator(y3[s==1],pik3[s==1],N,type="total")
cat("Case 4\n")
ht[4]=HTestimator(y4[s==1],pik4[s==1])
hajek[4]=Hajekestimator(y4[s==1],pik4[s==1],N,type="total")
ss[i,]=ht
ss1[i,]=hajek
}




###################################################
### code chunk number 6: up5
###################################################
#true values
tv=c(sum(y1),sum(y2),sum(y3),sum(y4))
for(i in 1:4)
{
cat("Case ",i,"\n")
cat("The mean of the Horvitz-Thompson estimators:",mean(ss[,i])," and the true value:",tv[i],"\n") 
MSE1=var(ss[,i])+(mean(ss[,i])-tv[i])^2
cat("MSE Horvitz-Thompson estimator:",MSE1,"\n")
cat("The mean of the Hajek estimators:",mean(ss1[,i])," and the true value:",tv[i],"\n") 
MSE2=var(ss1[,i])+(mean(ss1[,i])-tv[i])^2
cat("MSE Hajek estimator:",MSE2,"\n")
cat("Ratio of the two MSE:", MSE1/MSE2,"\n")
}



###################################################
### code chunk number 7: HT_Hajek_estimators.Snw:140-149 (eval = FALSE)
###################################################
## data(belgianmunicipalities)
## attach(belgianmunicipalities)
## # sample size
## n=20
## pik=inclusionprobabilities(Tot04,n)
## N=length(pik)
## 
## 
## sim=10
## ss=ss1=array(0,c(sim,4))
## 
## cat("Case 1\n")
## y1=rep(3,N)
## cat("Case 2\n")
## y2=TaxableIncome
## cat("Case 3\n")
## x=1:N
## pik3=inclusionprobabilities(x,n)
## y3=1/pik3
## cat("Case 4\n")
## epsilon=rnorm(N,0,sqrt(1/3))
## pik4=pik3
## y4=5*(x+epsilon)
## 
## ht=numeric(4)
## hajek=numeric(4)
## for(i in 1:sim)
## {
## cat("Simulation ",i,"\n")
## cat("Case 1\n")
## s=UPtille(pik)
## ht[1]=HTestimator(y1[s==1],pik[s==1])
## hajek[1]=Hajekestimator(y1[s==1],pik[s==1],N,type="total")
## cat("Case 2\n")
## s1=UPpoisson(pik)
## ht[2]=HTestimator(y2[s1==1],pik[s1==1])
## hajek[2]=Hajekestimator(y2[s1==1],pik[s1==1],N,type="total")
## cat("Case 3\n")
## ht[3]=HTestimator(y3[s==1],pik3[s==1])
## hajek[3]=Hajekestimator(y3[s==1],pik3[s==1],N,type="total")
## cat("Case 4\n")
## ht[4]=HTestimator(y4[s==1],pik4[s==1])
## hajek[4]=Hajekestimator(y4[s==1],pik4[s==1],N,type="total")
## ss[i,]=ht
## ss1[i,]=hajek
## }
## 
## 
## #true values
## tv=c(sum(y1),sum(y2),sum(y3),sum(y4))
## for(i in 1:4)
## {
## cat("Case ",i,"\n")
## cat("The mean of the Horvitz-Thompson estimators:",mean(ss[,i])," and the true value:",tv[i],"\n") 
## MSE1=var(ss[,i])+(mean(ss[,i])-tv[i])^2
## cat("MSE Horvitz-Thompson estimator:",MSE1,"\n")
## cat("The mean of the Hajek estimators:",mean(ss1[,i])," and the true value:",tv[i],"\n") 
## MSE2=var(ss1[,i])+(mean(ss1[,i])-tv[i])^2
## cat("MSE Hajek estimator:",MSE2,"\n")
## cat("Ratio of the two MSE:", MSE1/MSE2,"\n")
## }
## 
## 
## 
## sampling.newpage()
## 


