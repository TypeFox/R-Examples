bd.ME.optim.rho.all <-
function (x,t,sampling,yule = FALSE,maxitk=5,init=c(-1),posdiv=FALSE,survival=1)  {
	if (yule==TRUE) {
		print("not implemented. ext>0 !")
		} else {
	x<-sort(x)
    dev <- function(p) {
    	l<-vector()
    	mu<-vector()
    	help <- (length(p)+1)/3
    	for (i in 1:help) {
    		ltemp<-p[i+help]/(1-p[i])
    		mutemp<-p[i] * p[i+help]/(1-p[i])
			l<-c(l, ltemp)
    		mu<- c(mu,mutemp)
   	     }
    	rho<-c(sampling,p[(2*help+1):length(p)])
        out<- LikShifts(x,t,l,mu,rho,posdiv,survival) 
        out
    }
    numb<-length(t)
    if (init[1]==-1){
    init<-vector()
    for (i in 1:numb){
    	init<-c(0.001,init,0.05)
    	}
    init <- c(init,(2:numb*0+0.7))}
    out <- optim(init, dev,control=list(maxit=10000))#,method="SANN" )
    k<-1
    while (out$convergence != 0 && k<maxitk && 10000*10^k<.Machine$integer.max  ){ 
    		out <- optim(init, dev,control=list(maxit=10000*10^k)) 
    	 	k<-k+1
    	 }
	para<-1
	res<-list(out,para)
	}
	res
    }

