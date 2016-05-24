# This script performs the deviance calculations associated with Eq. 2 in our active MS
rm(list=ls(all=TRUE))
cols=c("city","sex","doseg","agexg","calg","kerma","PY","adjPY","num.entering",
"age","agex","tsx","cal","sv","gam","neut","lymphoma","NHL","leukemia","AML","ALL","CML","ATL","MM")      
d<-read.table("c:/data/abomb/HEMA87.dat", header=F,col.names=cols);
d=d[d$adjPY>0,] #remove two recs with zero py
d=d[d$kerma==1,] # take only kerma < 4 Gy
d=d[d$city==1,]#hiroshima=1, comment for pooling
#d=d[d$age>=20,]# restriction to adults => male difference regardless of city pooling
d$py=10^4*d$adjPY
d$calg=as.integer(cut(d$calg,c(0,2,4,6,8,10))) # pairwise binning of waiting time groups
m=d[d$sex==1,]; f=d[d$sex==2,] 
agem=55  # pooling cities holds significance with age>20 only if agem=55 centering is used
agem=0   # Hiroshima alone, with or without age>20, holds significance regardless of agem=0 or 55

require(bbmle)
nLL<-function(c0,k,beta,L1,L2,L3,L4,L5,x,agem) with(x, {L=c(L1,L2,L3,L4,L5)
				mn = (exp(c0+k*(age-agem)) + exp(-beta*abs(agex-30)/28.85)*sv*exp(L[calg]))*py
				-sum(stats::dpois(CML, mn, log=TRUE))})	
fit0 <- mle2(nLL,start=list(c0=-10,k=.04,L1=-10,L2=-10,L3=-10,L4=-10,L5=-10,beta=0.5),data=list(x=m,agem=agem))
fit0Fx <- mle2(nLL,start=list(c0=-10,k=.04,L1=-10,L2=-10,L3=-10,L4=-10,L5=-10),fixed = list(beta=0),data=list(x=m,agem=agem))
deviance(fit0Fx)-deviance(fit0)

fit0 <- mle2(nLL,start=list(c0=-10,k=.04,L1=-10,L2=-10,L3=-10,L4=-10,L5=-10,beta=0.5),data=list(x=f,agem=agem))
fit0Fx <- mle2(nLL,start=list(c0=-10,k=.04,L1=-10,L2=-10,L3=-10,L4=-10,L5=-10),fixed = list(beta=0),data=list(x=f,agem=agem))
deviance(fit0Fx)-deviance(fit0)

