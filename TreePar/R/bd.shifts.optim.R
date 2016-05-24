bd.shifts.optim <- function(x,sampling,grid,start,end,maxitk=5,yule=FALSE,ME=FALSE,all=FALSE,posdiv=FALSE,miniall=c(0),survival=1,groups=0){
	print("startest")
	print("test")
	x<-sort(x)	
	shifts<-length(sampling)-1
	estall<-list()
	convfail<-vector()
	if (length(miniall)==1){
		miniall<-list()
		est0<-bd.ME.optim(x,c(0),c(sampling[1]),yule,survival=survival,posdiv=posdiv,groups=groups)
		if (yule==FALSE){
		if (est0[[1]]$convergence != 0){
			print("convergence problem")
			convfail<-rbind(convfail,c(0,0))
		}}
		miniall<-c(miniall,list(c(est0[[1]]$value,est0[[1]]$par)))
		estall<-c(estall,list(est0)) 
		timeshifts<-c(0)
		lower<-1
		timevec<-vector()
	} else {
		lower<-length(miniall)
		help<-miniall[[lower]]
		timeshifts<-c(0,help[(length(help)-lower+2):length(help)])
	}
	if (length(sampling)>1){
			cuts<-round((end-start)/grid)	
	for (k in lower:shifts){
		minitime<-0
		est<-list()
		lik<-vector()
		timevec<-vector()
		mini<-c(200000)
		miniindex<-0
		jreal<-0
		for (j in 0:(cuts)){
			time1<-start + j*grid  # /cuts*(end-start)
			if (length(which(timeshifts==time1))==0){
			jreal<-jreal+1
			timevec<-c(timevec,time1)
			timetemp<-sort(c(timeshifts,time1))
			dupl<-order(c(timeshifts,time1))[length(c(timeshifts,time1))]-1
			if (ME==FALSE){
				inittemp<-miniall[[k]][2:((2*k)+1)]
				init<-c(inittemp[1:dupl],inittemp[dupl:(k+dupl)],inittemp[(k+dupl):length(inittemp)])
				}
			else if ( all==TRUE){
				if (k>1){
					inittemp<-miniall[[k]][2:((3*k))]
					if ((2*k+dupl)<=length(inittemp)){
						init<-c(inittemp[1:dupl],inittemp[dupl:(k+dupl)],inittemp[(k+dupl):(2*k+dupl-1)],1,inittemp[(2*k+dupl):length(inittemp)])
					}else{
						init<-c(inittemp[1:dupl],inittemp[dupl:(k+dupl)],inittemp[(k+dupl):(2*k+dupl-1)],1)}
				} else {
					inittemp<-miniall[[k]][2:((2*k)+1)]
					init<-c(inittemp[1:dupl],inittemp[dupl:(k+dupl)],inittemp[(k+dupl):length(inittemp)],1)
				}
			} else {
				inittemp<-miniall[[k]][2:(k+2)]
				if (length(inittemp)> (dupl+1)){
					init<- c(inittemp[1:(dupl+1)],1, inittemp[(dupl+2):length(inittemp)])
				} else {
					init<- c(inittemp[1:(dupl+1)],1)
				}}	
			if (ME==FALSE){
				temp<-bd.ME.optim(x,timetemp,sampling[1:length(timetemp)],yule,maxitk,init,posdiv,survival,groups=groups)}
				#test<-partransformvector(c(init,timetemp[-1]))
			else if ( all==TRUE){
				temp<-bd.ME.optim.rho.all(x,timetemp,sampling[1],yule,maxitk,init,posdiv,survival)
				print("higher taxa are ignored. groups = 0")
				#test<-partransformvectorrho(c(init,timetemp[-1]),sampling[1])
			} else {
				print(init)
				temp<-bd.ME.optim.rho(x,timetemp,sampling[1],yule,maxitk,init,posdiv,survival)
				print("higher taxa are ignored. groups = 0")}
			if (temp[[1]]$convergence != 0){
				print("convergence problem")
				convfail<-rbind(convfail,c(k,time1))
				}
			lik<-c(lik,temp[[1]]$value)
			if (temp[[1]]$value<mini[1] && length(which(time1==timeshifts))==0){
				mini<-c(temp[[1]]$value,temp[[1]]$par,timetemp[2:length(timetemp)])
				minitime<-timetemp
				miniindex<-(jreal+1)
			}
			#print(c("numbershifts","shift time","init"))
			print(c(k,time1))
			print(c(temp[[1]]$value,temp[[1]]$par,timetemp[2:length(timetemp)]))
			#print(test2)
			#print(init)
			est<-c(est,list(temp))	
		}
		}
	estall<-c(estall,list(est))
	miniall<-c(miniall,list(mini))
	print(miniall)
	timeshifts<-minitime
	}}
	out<-list(estall,miniall,timevec,convfail)
	out
}
