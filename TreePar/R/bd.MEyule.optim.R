bd.MEyule.optim <-
function (x,t,rho,maxitk=5,init=c(-1),groups=0)  {
	x<-sort(x)
    dev <- function(p) {
        out<- LikShifts(x,t,p,p*0,rho,groups)
        out
    }
    numb<-length(rho)
    if (init[1]== -1){
    init<-vector()
    for (i in 1:numb){
    	init<-c(init,1)
    	}}  else {
    init<-init[1:(length(init)/2)]		
    }
    if (length(init)==1){
	    out <- optimize(dev,interval=c(0,10))	
	    out$par<-out$minimum
	    out$value<-out$objective	

	} else {
	    	out <- optim(init, dev,control=list(maxit=10000))#,method="SANN"
	        k<-1
    		while (out$convergence != 0 && k<maxitk&& 10000*10^k<.Machine$integer.max ){ 
    			out <- optim(init, dev,control=list(maxit=10000*10^k)) 
    			k<-k+1
    	 	}
	}
	para<-1
	list(out,para)
    }

