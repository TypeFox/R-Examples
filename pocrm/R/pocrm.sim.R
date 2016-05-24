pocrm.sim<-function(r,alpha,prior.o,x0,stop,n,theta,nsim,tox.range){

sim <- sim1 <- apred <- lik <- pord <- ord <- ahat <- rpred <- next.lev <- n1 <- N <- NULL

d<-ncol(alpha)
s<-nrow(alpha)
q<-2

###LOAD FUNCTION 'crm'

crm<-function(obs,alpha,prior.o,theta){

sim<<-table(obs$level,obs$tox)
ifelse(dim(sim)[2]==1,sim1<<-data.frame(as.numeric(row.names(sim)),sim[,1],1-sim[,1]),sim1<<-data.frame(as.numeric(row.names(sim)),sim[,1],sim[,2]))
names(sim1)<<-c('level','nontox','tox')

apred<<-rep(0,s)
	lik<<-rep(0,s)
	for(k in 1:s){			
		ll<-function(a){
			la<-0 #value of the log-likelihood
			for(i in sim1$level){
			index<-match(i,sim1$level)
			la<-la+sim1$tox[index]*a*log(alpha[k,][i])+sim1$nontox[index]*log((1-alpha[k,][i]**a))
			}
			la
			}
		apred[k]<<-optimize(f=ll,interval=c(0,500),maximum=T)$maximum
		lik[k]<<-ll(apred[k])
		}
pord<<-(exp(lik)*prior.o)/sum(exp(lik)*prior.o)
#library("nnet") not necessary because listed as package dependency
ord<<-which.is.max(pord)
ahat<<-apred[ord]
rpred<<-alpha[ord,]**ahat
next.lev<<-which.is.max(-(abs(rpred-theta)))
next.lev
    }
###'crm' ENDS HERE


###LOAD FUNCTION 'twostgcrm'

twostgcrm<-function(r,x0,stop,n,theta){

	n1<<-n+1
	obs<-data.frame(cbind(1:n1,rep(0,n1),rep(0,n1),rep(0,n1),rep(0,n1)))
	names(obs)<-c('patient','level','tox','a','order')

#######Beginning at dose level 1
###1st Stage:  Up and down scheme with cohort size 1

i<-1
#x0<-lapply(zones,ff)
##'initial.scheme' is a vector indicating the Stage I escalation scheme
initial.scheme<-x0
initial.scheme<-c(initial.scheme,rep(ncol(alpha),n1-length(initial.scheme)))

while(i < n1){
obs$order[i]<-99
obs$level[1]<-initial.scheme[1]
p<-runif(1)
#number of tox in 1st patient
index<-p<=r[obs$level[i]] ##determines any toxicities
obs$tox[i]<-obs$tox[i]+as.numeric(index)
if(any(obs$tox[1:i]==1) & any(obs$tox[1:i]==0)){
	q<-2
	break
}
if(all(obs$tox[1:i]==1)){
	i<-i+1
	obs$level[i]<-initial.scheme[1]
}
if(all(obs$tox[1:i]==0)){
	i<-i+1
	obs$level[i]<-initial.scheme[i]
}

if(length(obs$level[obs$level==d])==stop+1){
	MTD<-d
	break
}
}
	

##2nd stage
N<<-table(obs$level>0)[2]+1
if(any(obs$tox>0)){
		  level<-crm(obs[1:(N-1),],alpha,prior.o,theta)
		  obs$a[N-1]<-ahat
	        obs$order[N-1]<-ord

for(j in N:n1){
##assigment for remaining patients
	obs$level[j]<-level
	if(obs$level[n1]>0){
			MTD<-obs$level[n1]
		break
		}

	if(length(obs$level[obs$level==level])==stop+1){
		MTD<-level
		break
		}

	index<-runif(1)<=r[obs$level[j]]
	if(index){obs$tox[j]<-1}
	level<-crm(obs[1:j,],alpha,prior.o,theta)
	obs$a[j]<-ahat
	obs$order[j]<-ord
	
##crm dose recommendation for Nth patient
	}
	} else
	MTD<-d
	out<-list(trial=obs[obs$level>0,],MTD.selection=MTD)
	}
###'twostgcrm' ENDS HERE


##LOAD FUNCTION 'pocrm.simul'

pocrm.simul<-function(nsim,tox.range){
dlt<-p<-ca<-rep(0,d)
rec0<-ss<-rep(0,nsim)
a<-b<-matrix(,nrow=d,ncol=nsim)
for(g in 1:nsim){
	fit<-twostgcrm(r,x0,stop,n,theta)
	ss[g]<-sum(fit$trial$order>0)
	   for(w in 1:d){
			a[w,g]<-sum(fit$trial$level[1:ss[g]]==w)
			b[w,g]<-sum(fit$trial$level[1:ss[g]]==w & fit$trial$tox[1:ss[g]]==1)		
			rec0[g]<-fit$MTD.selection
			}
}
for(h in 1:d){
	p[h]<-sum(rec0==h)/nsim
	ca[h]<-sum(a[h,])/sum(ss)
}
output<-list(true.prob=r,MTD.selection=round(p,2),patient.allocation=round(ca,2),percent.DLT=sum(b)/sum(ss),mean.n=mean(ss),acceptable=sum(p[which(round(abs(r-theta),2)<=tox.range)]))
}
##'pocrm.sim' END HERE

if(nsim==1){
	out<-twostgcrm(r,x0,stop,n,theta)
} else{
	out<-pocrm.simul(nsim,tox.range)
}
out
}

