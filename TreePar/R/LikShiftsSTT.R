#numbd: either as many as birth rates (default:numbd=0); or just numbd=1.
#if tconst vector: fixed shift times; if -1 not fixed
#sprob: probability of sampling an extinct individual
LikShiftsSTT <- function(par,times,ttype,numbd=0,tconst=-1,sampling=0,sprob,root=0,survival=1,tfixed=vector(),mint=0,maxt=0) {
	if (root == 0 && sum(ttype)*2 != length(ttype)) {print("Invalid input. Most likely you did not add a root edge to your tree using addroot().")} else {
	rho<-sampling	
	if (length(ttype)==1 && ttype==0){ttype<-times*0+1}
	shifts<-length(sprob)
	if (maxt==0) {maxt<-max(times)}	
	
	#set parameters
	if (shifts>1){
		l<- par[1:shifts] 
		k<-shifts+1
		if (numbd==1) {mu<-(1:shifts)*0+par[k]*(1-sprob)
			psi<-(1:shifts)*0+par[k]*(sprob)
			k<-k+1
		} else {
			mu<-par[k:(k+shifts-1)]*(1-sprob)
			psi<-par[k:(k+shifts-1)]*(sprob)
			k<-k+shifts
		}
		if (tconst == -1){t<-sort(c(0,par[k:length(par)],tfixed))} else {t<-c(0,tconst)}
	} else {
		rho<-c(rho,0)
		t<-c(0,0)
		mu<-c(par[2],par[2])*(1-sprob)
		psi<-c(par[2],par[2])*sprob
		l<-c(par[1],par[1])
	}
	
	if (length(rho)==1 && rho == 0){rho<-l*0}
		
	transmission<-times[which(ttype==1)]
	sampling<-times[which(ttype==0)]  
	samplingkeep<-sampling
	if (root==1){
		transmission<-c(transmission, max(transmission))
	}
	extant<-length(transmission)-length(sampling)
	out<- -10^12
	boundary<-0
	#check if valid parameters
	for (i in 1:length(l)){	if (l[i]<=0 || min(psi[i],mu[i])<0 ){boundary<-1}}
	#if (shifts==2 && maxt==max(times)) {if (t[2]<sort(transmission)[5] || t[2]>sort(transmission,decreasing=TRUE)[5])  boundary<-1}	
	if (shifts==2 && maxt<max(times)) {if (t[2]<=mint || t[2]>=maxt)  boundary<-1}	
	if (shifts>2) {
			if (par[length(par)]<=mint || par[length(par)]>=maxt)  boundary<-1
	}		

	if (boundary==0) {
		out<- -(root+1)*log(2*l[interstt(max(transmission),t)])
		if (survival==1){
			index<- interstt(max(transmission),t)
			out<- out - (root+1)* log(1- p(index,max(transmission),t,l,mu,psi,rho))
			#print("1-Ext")
			#print(1- p(index,max(transmission),t,l,mu,psi,rho))
		}

		if (extant>0){out<- out+ extant*log(rho[1]) }
		#all the samples at time t
		if (t[2]!=0) {
			for (j in 2:length(t)){
				index<-which(sampling==t[j])
				if (length(index)>0){
				sampling<-sampling[-index]
				out<-out+length(index)*log(rho[j])}
			}}

	
		for (j in 1:length(transmission)){
			out<-out+(qfuncskylog(transmission[j],t,l,mu,psi,rho))+log(2*l[interstt(transmission[j],t)])
			#print("q at t")
			#print(c(6-transmission[j],exp(qfuncskylog(transmission[j],t,l,mu,psi,rho))))
		}
	#sampling<-times[which(ttype==0)] #check


		if (length(sampling)>0){
			for (j in 1:length(sampling)){
				out<-out-(qfuncskylog(sampling[j],t,l,mu,psi,rho))+log(psi[interstt(sampling[j],t)])
				#print("q at t")
				#print(c(6-sampling[j],exp(qfuncskylog(sampling[j],t,l,mu,psi,rho))))

		}}
		for (j in 2:length(t)) {
			samplingtemp<-extant
			if (length(samplingkeep)>0){
				samplingtemp<-samplingtemp+length(which(samplingkeep<t[j]))}
			transmissiontemp<-length(which(transmission<t[j]))
			nj<-samplingtemp-transmissiontemp
			out<-out+nj*(qfuncskylog(t[j],t,l,mu,psi,rho))
			#print("q at t")
			#print(c(6-t[j],exp(qfuncskylog(t[j],t,l,mu,psi,rho)),nj))

		}			
	}
	out<- out-(length(transmission)-1-root)*log(2)
	-out
	}
}
