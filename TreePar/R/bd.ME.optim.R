bd.ME.optim <-
function (x,t,rho,yule,maxitk=5,init=c(-1),posdiv=FALSE,survival=1,groups=0)  {
	if (yule==TRUE) {
		res<-bd.MEyule.optim(x,t,rho,groups)
		} else {
	x<-sort(x)
    dev <- function(p) {
    	boundary=1
    	l<-vector()
    	mu<-vector()
    	help <- length(p)/2
    	for (i in 1:help) {
    		ltemp<-p[i+help]/(1-p[i])
    		mutemp<-p[i] * p[i+help]/(1-p[i])
			l<-c(l, ltemp)
    		mu<- c(mu,mutemp)
    	}
      	out<- LikShifts(x,t,l,mu,rho,posdiv,survival,groups) 
        out
    }
    numb<-length(rho)
    if(init[1] == -1){
    init<-vector()
    for (i in 1:numb){
    	init<-c(0.001,init,0.05)
    }}
	out <- try(optim(init, dev,control=list(maxit=10000)),silent=TRUE) #,method="SANN"
	if (class(out) == "try-error"){init<-10*init
	out <- try(optim(init, dev,control=list(maxit=10000)),silent=TRUE)} #,method="SANN"
	if (class(out) == "try-error"){init<-init/100
	out <- optim(init, dev,control=list(maxit=10000))} #,method="SANN"
	k<-1
    while (out$convergence != 0 && k<maxitk && 10000*10^k<.Machine$integer.max ){ 
    	out <- optim(init, dev,control=list(maxit=10000*10^k)) 
    	k<-k+1
    }
	para<-1
	res<-list(out,para)
	}
	res
    }

