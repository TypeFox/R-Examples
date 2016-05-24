pocrm.imp <-
function(alpha,prior.o,theta,y,combos){
	
	data<-as.matrix(table(combos,y))
	level<-as.numeric(row.names(data))
	nontox<-as.numeric(data[,1])
	tox<-as.numeric(data[,2])

apred<-rep(0,nrow(alpha))
	lik<-rep(0,nrow(alpha))
	for(k in 1:nrow(alpha)){			
		ll<-function(a){
			la<-0 #value of the log-likelihood
			for(i in level){
			index<-match(i,level)
			la<-la+tox[index]*a*log(alpha[k,][i])+nontox[index]*log((1-alpha[k,][i]**a))
			}
			la
			}
		apred[k]<-optimize(f=ll,interval=c(0,500),maximum=T)$maximum
		lik[k]<-ll(apred[k])
		}
pord<-(exp(lik)*prior.o)/sum(exp(lik)*prior.o)
#library("nnet") not necessary because listed as package dependency
ord<-which.is.max(pord)
ahat<-apred[ord]
rpred<-alpha[ord,]**ahat
next.lev<-which.is.max(-(abs(rpred-theta)))
out<-list(ord.prob=round(pord,3),order.est=ord,a.est=round(ahat,3),ptox.est=round(rpred,3),dose.rec=next.lev)
}
