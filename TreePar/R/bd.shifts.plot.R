bd.shifts.plot<-function(resall,shifts,timemax=100,ratemin=-1,ratemax=1,plotturnover=FALSE) {
	pick<-function(resall,shifts){
		resall[[shifts+1]]
	}
	estimatesall<-sapply(resall,pick,shifts=shifts)
	plot(c(0,-timemax),c(ratemin,ratemax),col="white",xlab="time before the present",ylab="diversification rate")
	
	#PROBLEM IF ONLY ONE SHIFT!!!
	
	for (i in 1:length(estimatesall[1,])){
		estimates<-estimatesall[,i]
		rates<-length(estimates)/3
		estimates<-estimates[-1]
		if (rates>1){
			time<-estimates[(length(estimates)-rates+2):length(estimates)]
			time<-sort(c(time,time,0,timemax))
			turnover<-estimates[1]
			div<-estimates[rates+1]
			for (j in 1:(rates-1)){turnover<-c(turnover,estimates[j:(j+1)])
				div<-c(div,estimates[(rates+j):(rates+j+1)])
			}
			turnover<-c(turnover,estimates[rates])
			div<-c(div,estimates[2*rates])}
		else {time<-c(0,timemax)
			turnover<-c(estimates[1],estimates[1])
			div<-c(estimates[2],estimates[2])
		}
	if (plotturnover==TRUE){lines(-time,turnover,col="red")}
	lines(-time,div,col="blue")
	}
}

