# confInts.R  
# Exponential background CML parameter estimate confidence intervals

# This script starts out by using optim() and an explicitly defined
# log-likelihood function to fit exp(c0+k*age) to CML incidence for ages >20 
# using Poisson regression.  It then shows that the same results are obtained by
# either mle2() or by glm(), but slightly different results are obtained using 
# the least squares functions nls() and lm().  This script also shows that when 
# ages are centered about a median age of 55, the confidence interval of c0 (the
# log-space intercept) is shortened more so than that of k (the log-space 
# slope).

# NOTICE: you must
if (0) {  # switch this to 1 (i.e. run this chunk) if you never ran it before
	library(SEERaBomb)
	(df=getFields())
	(df=pickFields(df))
	mkSEERold(df,dataset="00")  
}


rm(list=ls(all=TRUE))
seerHome="~/data/SEER"
load(file.path(seerHome,"00/pops.RData")) # this loads in pops
pym=NULL;pyf=NULL
for (i in 0:18) 
{   pym[i+1]=with(pops,sum(population[(popsex==1)&(popage==i)]))
	pyf[i+1]=with(pops,sum(population[(popsex==2)&(popage==i)])) }

load(file.path(seerHome,"00/lymyleuk.RData")) # this loads in DF
DF$agerec=as.numeric(cut(DF$agedx,breaks=c(0,1,seq(5,85,5),130)))
head(DF)
d=DF[(DF$cancer=="CML")&(DF$seqnum<=1),] # paper used SEER 1973-2008
m=hist(d$agerec[d$sex==1],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
f=hist(d$agerec[d$sex==2],breaks=c(seq(-.5,17.5,1),100),plot=FALSE)$counts
age=c(0.5,3,seq(7.5,87.5,5))
agem=55  # center ages in fits
agem=0
datam=data.frame(age=age-agem,cases=m,py=pym)[6:19,] 
dataf=data.frame(age=age-agem,cases=f,py=pyf)[6:19,]

fexp<-function(p,dat)	{	c0=p[1];k=p[2]; 
	with(dat,{mn=exp(c0+k*age)*py
				-sum(cases*log(mn) - mn)})	}
p0=c(c0=-8,k=.04)
ssolm=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=datam,hessian=TRUE)
sig=sqrt(diag(solve(ssolm$hessian)))
(CIM=cbind(point=ssolm$par,lower=ssolm$par-1.96*sig,upper=ssolm$par+1.96*sig))
CIM[,3]-CIM[,2] # The intercept is ~3-fold more accurate when the data is centered around it!
p0=c(c0=-10,k=.04)
ssolf=optim(p0,fexp,method="L-BFGS-B",control=list(maxit=400),dat=dataf,hessian=TRUE)
sig=sqrt(diag(solve(ssolf$hessian)))
(CIF=cbind(point=ssolf$par,lower=ssolf$par-1.96*sig,upper=ssolf$par+1.96*sig))

# now compare CI by different fitting methods
# mle2 is an optim wrapper in bbmle that has some conveniences 
require(bbmle) # note that dataframes used globally below are already centered in age
nLL<-function(c0,k,x) with(x, 	
	{mn=exp(c0+k*age)*py
	-sum(stats::dpois(cases, mn, log=TRUE))})
fit0m <- mle2(nLL,start=list(c0=-10,k=.04),data=list(x=datam))
fit0f <- mle2(nLL,start=list(c0=-10,k=.04),data=list(x=dataf))
glm1m=glm(cases~age+offset(log(py)),data=datam,family=poisson)
glm1f=glm(cases~age+offset(log(py)),data=dataf,family=poisson)
datam=transform(datam,incid=cases/py)
dataf=transform(dataf,incid=cases/py)
nl1m=nls(incid~exp(c0+k*age),data=datam,start=c(c0=-10,k=0.1))
nl1f=nls(incid~exp(c0+k*age),data=dataf,start=c(c0=-10,k=0.1))
lm1m=lm(log(incid)~age,data=datam)
lm1f=lm(log(incid)~age,data=dataf)

CIM
cbind(coef(glm1m),confint(glm1m))
cbind(coef(fit0m),confint(fit0m))
cbind(coef(lm1m),confint(lm1m))
cbind(coef(nl1m),confint(nl1m))

CIF
cbind(coef(glm1f),confint(glm1f))
cbind(coef(fit0f),confint(fit0f))
cbind(coef(lm1f),confint(lm1f))
cbind(coef(nl1f),confint(nl1f))

# conclusions: optim, mle2 and glm all give the same result
# since they are all Poisson regression. The least squares
# approaches can yield noticeably different results.

