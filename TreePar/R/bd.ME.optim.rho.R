bd.ME.optim.rho <-
function (x,t,sampling,yule = FALSE,maxitk=5,init=c(-1),posdiv=FALSE,survival=1)  {
	if (yule==TRUE) {
		print("not implemented. ext>0 !")
		} else {
	x<-sort(x)
    dev <- function(p) {
    	rho<-c(sampling,p[3:length(p)])
    	l<-(1:length(rho))*0+ (p[2]/(1-p[1]))
    	mu<-(1:length(rho))*0+  (p[1] * p[2]/(1-p[1]))
      	out<- LikShifts(x,t,l,mu,rho,posdiv,survival)
        out
    }
    numb<-length(t)
    if (init[1]==-1){
    init<-vector()
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

