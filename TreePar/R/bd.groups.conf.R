bd.groups.conf <- function(out,dev)    {
	para <- out$par
    #se<-c(NA,NA)
    #try(se <- sqrt(diag(solve(out$hessian))))
    k<- length(out$par)
    Dev <- out$value
    CI <- matrix(NA, k, 2)
    
    for (i in 1:k) {
    	if (i==1) {
    		foo <- function(p) dev(c(p, para[(i+1):k])) - 3.84/2 - Dev
    	} else if (i==k) {
    		foo <- function(p) dev(c(para[1:(k-1)],p)) - 3.84/2 - Dev
    	} else {
    		foo <- function(p) dev(c(para[1:(i-1)],p, para[(i+1):k])) - 3.84/2 - Dev
    	}
    	inc <- 0.01
    	lo <- para[i] - inc
    	up <- para[i] + inc
    	while (foo(lo) < 0) lo <- lo - inc
    	while (foo(up) < 0) up <- up + inc
    	CI[i, 1] <- uniroot(foo, lower = lo, upper = para[i])$root
    	if (CI[i, 1] < 0) 
        	CI[i, 1] <- 0
    	CI[i, 2] <- uniroot(foo, lower = para[i], upper = up)$root
    }
    r<- vector()
    d<-vector()
    t<-vector()
    for (i in 1:((k+1)/3)) {
    	r<-c(r,paste("d",i-1,"/b",i-1,sep=""))
    	d<-c(d,paste("b",i-1,"-d",i-1,sep=""))
    	t<-c(t,paste("t",i))
    	}
    names <- c(r,d,t)
    names<-names[-length(names)]
    CI<-cbind(out$par,CI)#,se)
    rownames(CI)<-names
    colnames(CI) <- c("MLE","lo", "up")#,"se")
    CI
	}