### R code from vignette source 'UPexamples.Snw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: UPexamples.Snw:21-25
###################################################
library(sampling)
ps.options(pointsize=12)
options(width=60)



###################################################
### code chunk number 2: entropy1
###################################################
data(belgianmunicipalities)
attach(belgianmunicipalities)
n=50


###################################################
### code chunk number 3: entropy2
###################################################
pik=inclusionprobabilities(averageincome,n)


###################################################
### code chunk number 4: entropy3
###################################################
s=UPmaxentropy(pik)


###################################################
### code chunk number 5: entropy4
###################################################
as.character(Commune[s==1])


###################################################
### code chunk number 6: entropy5
###################################################
pi2=UPmaxentropypi2(pik)


###################################################
### code chunk number 7: entropy6
###################################################
rowSums(pi2)/pik/n


###################################################
### code chunk number 8: entropy7
###################################################
data(belgianmunicipalities)
attach(belgianmunicipalities)
pik=inclusionprobabilities(averageincome,50)
pik=pik[pik!=1]
n=sum(pik)
pikt=UPMEpiktildefrompik(pik)
w=pikt/(1-pikt)
q=UPMEqfromw(w,n)


###################################################
### code chunk number 9: entropy8
###################################################
UPMEsfromq(q)



###################################################
### code chunk number 10: entropy9
###################################################
sim=10000
N=length(pik)
tt=rep(0,N)
for(i in 1:sim) tt = tt+UPMEsfromq(q)
tt=tt/sim
max(abs(tt-pik))


###################################################
### code chunk number 11: up1
###################################################
b=data(belgianmunicipalities)
pik=inclusionprobabilities(belgianmunicipalities$Tot04,200)
N=length(pik)
n=sum(pik)


###################################################
### code chunk number 12: up2
###################################################
sim=10
ss=array(0,c(sim,8))


###################################################
### code chunk number 13: up3
###################################################
y=belgianmunicipalities$TaxableIncome


###################################################
### code chunk number 14: up4
###################################################
ht=numeric(8)
for(i in 1:sim)
{
cat("Step ",i,"\n")
s=UPpoisson(pik)
ht[1]=HTestimator(y[s==1],pik[s==1])
s=UPrandomsystematic(pik)
ht[2]=HTestimator(y[s==1],pik[s==1])
s=UPrandompivotal(pik)
ht[3]=HTestimator(y[s==1],pik[s==1])
s=UPtille(pik)
ht[4]=HTestimator(y[s==1],pik[s==1])
s=UPmidzuno(pik)
ht[5]=HTestimator(y[s==1],pik[s==1])
s=UPsystematic(pik)
ht[6]=HTestimator(y[s==1],pik[s==1])
s=UPpivotal(pik)
ht[7]=HTestimator(y[s==1],pik[s==1])
s=srswor(n,N)
ht[8]=HTestimator(y[s==1],rep(n/N,n))
ss[i,]=ht
}


###################################################
### code chunk number 15: up5
###################################################
colnames(ss) <- 
c("poisson","rsyst","rpivotal","tille","midzuno","syst","pivotal","srswor")
boxplot(data.frame(ss), las=3)



###################################################
### code chunk number 16: UPexamples.Snw:162-170 (eval = FALSE)
###################################################
## b=data(belgianmunicipalities)
## pik=inclusionprobabilities(belgianmunicipalities$Tot04,200)
## N=length(pik)
## n=sum(pik)
## sim=10
## ss=array(0,c(sim,8))
## y=belgianmunicipalities$TaxableIncome
## ht=numeric(8)
## for(i in 1:sim)
## {
## cat("Step ",i,"\n")
## s=UPpoisson(pik)
## ht[1]=HTestimator(y[s==1],pik[s==1])
## s=UPrandomsystematic(pik)
## ht[2]=HTestimator(y[s==1],pik[s==1])
## s=UPrandompivotal(pik)
## ht[3]=HTestimator(y[s==1],pik[s==1])
## s=UPtille(pik)
## ht[4]=HTestimator(y[s==1],pik[s==1])
## s=UPmidzuno(pik)
## ht[5]=HTestimator(y[s==1],pik[s==1])
## s=UPsystematic(pik)
## ht[6]=HTestimator(y[s==1],pik[s==1])
## s=UPpivotal(pik)
## ht[7]=HTestimator(y[s==1],pik[s==1])
## s=srswor(n,N)
## ht[8]=HTestimator(y[s==1],rep(n/N,n))
## ss[i,]=ht
## }
## colnames(ss) <- 
## c("poisson","rsyst","rpivotal","tille","midzuno","syst","pivotal","srswor")
## boxplot(data.frame(ss), las=3)
## 
## 
## sampling.newpage()
## 


