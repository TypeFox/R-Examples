HWEmodelsMLE <-
function(nvec){
       k <- .5*(-1+sqrt(1+8*length(nvec)))
       n <- matrix(0,nrow=k,ncol=k)
       count <- 1
       for(i in 1:k){
       	     for(j in 1:k){
                 ## n is lower triangular
                 if(i <= j){ n[i,j] <- nvec[count]; count<-count+1}
	     }
       }
       mu <- n; mt <- t(n)

###n2 <- matrix(0, nrow=k, ncol=k)
###n2[lower.tri(n2, diag=TRUE)] <- nvec

       ##      m is the number of alleles matrix
       m <- mt+mu; qhat <- fqhat <- NULL
       phat <- fixindhat <- matrix(0,nrow=k,ncol=k)
       for (i in 1:k) {
       	   qhat[i] <- sum(m[i,])/(2*sum(n))
       	   for (j in 1:k){
               ifelse(i==j,phat[i,j] <- .5*m[i,j]/sum(n),
                           phat[i,j] <- m[i,j]/sum(n))
           }
	}
	for(i in 1:k){
              cat("Allele Marginal prob: ",i,qhat[i],"\n")
              for(j in 1:k){
	      	    if (i>j){
		       fixindhat[i,j] <- 1-phat[i,j]/(2*qhat[i]*qhat[j])}}}
        init <- baselogit(rep(1/k,k))$baselogit
	maximum <- optim(par=c(init,0),fn=MultLogLik,nvec=nvec,control=list(fnscale=-1))
	cat("Convergence = ",maximum$conv,"(0 is successful convergence)\n")
	fmaxloglik <- maximum$value
	fqhat <- invbaselogit(maximum$par[1:k-1])$probs
	fsingle <- (exp(maximum$par[k])-1)/(exp(maximum$par[k])+1)
#	cat("Maximum LogLik (single f model) = ",fmaxloglik,"\n")
	pmin <- min(fqhat)
	fmin <- -pmin/(1-pmin)
	cat("Probs and f at max and fmin:\n ")
	cat(fqhat,fsingle,fmin,"\n")
	list(phat=phat,qhat=qhat,fmaxloglik=fmaxloglik,fqhat=fqhat,fsingle=fsingle,fmin=fmin,fixindhat=fixindhat)
}

