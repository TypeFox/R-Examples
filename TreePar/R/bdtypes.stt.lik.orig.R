#Calculate the likelihood numerically
#par must be  lambda11 lambda12 lambda21 lambda22 death1 death2 gamma12 gamma21
#states[i] belongs to leaf i
#neg value -m in fix mean we equal to the parameter m times the entry in thrid row
bdtypes.stt.lik.orig <- function(par,phylo,fix=rbind(c(0,0),c(0,0)),sampfrac,survival=0,posR=0,unknownStates=FALSE,root=0){  
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
		lik<-try(BDSSnum.help(phylo,1,l,m,psi,summary,unknownStates))
		if (class(lik)!="try-error"){
		LambMu<-l[1]-l[4]-(m[1]+psi[1])+(m[2]+psi[2])
		c<- sqrt(LambMu^2 +4*l[2]*l[3])
		f1<- (c+LambMu)/(c+LambMu+2*l[2])
		out<- try(-log((lik[3]*f1)/(1-lik[1])^(survival) + (lik[4]*(1-f1))/(1-lik[2])^(survival)))
		if((class(out)!="numeric") || (out=="NaN") ||  (out=="Inf") ){out <- 10^10}}
	}
	if (out==10^10){out<-10^1000}		#Tanja 24.4.2012
	if (prpar==TRUE){print(par)}
	out	
}
