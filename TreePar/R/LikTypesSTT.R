#Calculate the likelihood numerically
#par must be  lambda11 lambda12 lambda21 lambda22 death1 death2 gamma12 gamma21
#states[i] belongs to leaf i
#neg value -m in fix mean we equal to the parameter m times the entry in thrid row
LikTypesSTT <- function(par,phylo,fix=rbind(c(0,0),c(0,0)),sampfrac,survival=0,posR=0,unknownStates=FALSE,rtol=1e-12,atol=1e-12,migr=0,freq=0,cutoff=10^12){  
	root<-0 #root =1 not tested
	prpar<-FALSE
	maxpar<-100
	partemp<-vector()
	k<-1
	for (i in 1:8){
		index<-which(i == fix[1,])
		if (length(index)>0) {
			if (fix[2,index]>=0){
			partemp<-c(partemp,fix[2,index])} else {
			temp<- - fix[2,index]
			if (temp == 0.4){		#make lambdas in same ratio
				partemp<-c(partemp,partemp[3]*partemp[2]/partemp[1])	
				  }else {
			partemp<-c(partemp,partemp[temp]*fix[3,index])}
			}
		} else {
			partemp<-c(partemp,par[k])
			k<-k+1
		}
	}
	#print(partemp)
	
	death<-partemp[5:6]
	l<-partemp[1:4]
	gamma<-partemp[7:8]
	psi<-death*sampfrac
	m<-death*(1-sampfrac)	
	
	if (root==1){
		cut<-phylo$edge[1,1]
		for (i in 1:length(phylo$edge[,1])){
			if (phylo$edge[i,1] >= cut){phylo$edge[i,1]<-phylo$edge[i,1]+1}
			if (phylo$edge[i,2] >= cut){phylo$edge[i,2]<-phylo$edge[i,2]+1}
		}
		phylo$edge<-rbind(c(cut,phylo$edge[1,1]),phylo$edge)
		phylo$edge.length<-c(0,phylo$edge.length)
		}

	summary<-get.times2(phylo)
	out <-  10^10
	temp<-1
	R0temp<-try(R0types(l[1],l[2],l[3],l[4],death[1],death[2]))
	if (posR==1 && class(R0temp)=="numeric" && R0temp<1) {temp<-0}
	if (posR==1 && class(R0temp)=="try-error") {temp<-0}
	check<-((length(which(partemp=="NaN"))>0)||(min(l,psi))<0 || m<0  || max(l,m,psi)>maxpar || (temp==0))  #19.4.12: (min(l,psi))<=0
	if (check){out <-  10^10} else {   
		lik<-try(BDSSnum.help(phylo,1,l,m,psi,summary,unknownStates,rtol,atol,migr,cutoff))
	##	print("lik")
	##print(lik)
		if (class(lik)!="try-error"){
		LambMu<-l[1]-l[4]-(m[1]+psi[1])+(m[2]+psi[2])
		c<- sqrt(LambMu^2 +4*l[2]*l[3])
		f1<- (c+LambMu)/(c+LambMu+2*l[2])
		if (migr==1){
			LambMu<-l[1]-l[4]-(m[1]+l[2]+psi[1])+(m[2]+psi[2]+l[3]) #Phil Trans paper: l same, m->m+gamma
			c<- sqrt(LambMu^2 +4*l[2]*l[3])
			f1<- (c+LambMu)/(c+LambMu+2*l[2])
			}
		if (freq>0) {f1<-freq}

		#print(log((lik[3]*f1)/(1-lik[1])^(survival)))
		#print(lik[4]*(1-f1))
		#print(1-lik[2])
		out1<-0
		out2<-0
		
		# if ((1-lik[1])>10^(-40)){
			# out1<-try((lik[3]*(f1))/(1-lik[1])^(survival))}
		# if ((1-lik[2])>10^(-40)){out2<-try((lik[4]*(1-f1))/(1-lik[2])^(survival))}
		# out<- try(-log(out1 + out2))
		out <- ( lik[3]*(f1) + lik[4]*(1-f1)  )/(1-(f1*lik[1]+(1-f1)*lik[2]))^(survival) #Tanja July 2014 P[T|Surv]=P[T,S]/P[S]=(P[T1]+P[T2])/P[S]
		if((class(out)!="numeric") || (out=="NaN") ||  (out=="Inf") ){out <- 10^10}
		}
	}
	if (out==10^10){out<-10^1000}		#Tanja 24.4.2012
	if (prpar==TRUE){print(par)}
	##print("lik")
	##print(lik)
	-log(out)	
}
