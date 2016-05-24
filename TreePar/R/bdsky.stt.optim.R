#numbd=1: constant death rate the whole time
#else: varying death rate
#root=1: tree stops at mrca. root=0: also root edge.
bdsky.stt.optim <- function(x,ttype=0,rho=0,sampprob=c(0),constdeath=0,root=0) {
	numbd<-constdeath
	if (length(ttype)==1 && ttype==0) {ttype<-x*0+1}		#everything is sampling
	sprob<-sampprob
	times<-x
	shifts<-length(sprob)
	resminall<-vector()
	resall<-vector()
	likall<-vector()
	convissue<-vector()
	restemp <- subplex(c(2,1),BDSSanal,times=times,ttype=ttype,rho=rho,sprob=sprob[1],root=root)
	initrates<-restemp$par
	initt<-vector()
	likmin<-restemp$value
	#resmin<-restemp  #Tanja added July 16 2014
	#tempmin<-1		#Tanja added July 16 2014
	resminall<-c(resminall,list(c(restemp$value,restemp$par[2]/restemp$par[1],restemp$par[1]-restemp$par[2])))
	mininterval<-vector()
	if (shifts>1) {for (i in 2:shifts){		
		timesall<-vector()
		timessort<-sort(unique(times))
		intervals<-length(timessort)-1
		for (k in 2:intervals){
			if ((length(which(mininterval==k)))==0) {
			start<-timessort[k-1]
			end<-timessort[k]
			if (start<end){
				inittime<-start+(end-start)/2
				repl<-order(c(inittime,initt))[1]
				if (numbd==0){
					init<-c(initrates[1:repl],initrates[repl:(repl+length(initrates)/2)],initrates[(repl+length(initrates)/2):length(initrates)])} else { 
					init<-c(initrates[1:repl],initrates[repl:length(initrates)])} 
				restemp <- try(optim(c(init,inittime), LikShiftsSTT,times=times,ttype=ttype,numbd=numbd,sampling=rho,sprob=sprob[1:i],tfixed=initt,mint=start,maxt=end,root=root,control=list(maxit=10000),method="BFGS"))
				
				#method="BFGS",control=list(maxit=100000,reltol=10^(-15),abstol=10^(-15))))#,method="BFGS"))
				                   #,method="SANN")				#,control=list(tmax=50,temp=50,maxit=50000))
				print(paste("calculated", round(100*k/intervals,1), "% for",i,"rates"))
				if (class(restemp) != "try-error"){
				#print(c(restemp$value,restemp$par,restemp$convergence,inittime))
				if (likmin>restemp$value) {
					likmin<-restemp$value
					resmin<-restemp
					tempmin<-k
				}
				resall<-c(resall,list(restemp))
				likall<-rbind(likall,c(restemp$value,restemp$par[length(restemp$par)]))
				if (restemp$convergence != 0){convissue<-rbind(convissue,c(restemp$convergence,shifts-1,k))}} else {convissue<-rbind(-100,c(1,shifts-1,k))}
		}
		}}
		ltemp<-resmin$par[1:i]
		mutemp<-resmin$par[(i+1):(length(resmin$par)-1)]
		#print(ltemp)
		#print(mutemp)
		mininterval<-c(mininterval,tempmin)
		timetemp<-resmin$par[length(resmin$par)]
		tempmin<-c(resmin$value,mutemp/ltemp,ltemp-mutemp,sort(c(initt,timetemp)))
		resminall<-c(resminall,list(tempmin))
		print("resmin")
		print(resmin)
		initrates<-resmin$par[1:((length(resmin$par)-1))]
		initt<-c(initt,resmin$par[length(resmin$par)])
		#print(initrates)
		#print(initt)
		}}
	out<-list(resminall,convissue)
	out
}
