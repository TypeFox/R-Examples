gradient<-function(theta, w, x, beta, y, k, levels, Ymat, link){
    alpha <- theta[1:(k-1)]
	if (dim(w)[2]!=0) {
    	zeta <- theta[k:length(theta)]
		if (is.null(x)) {
    		Xb <- w%*%zeta
    	} else if (!is.null(x)) {
    		Xb <- w%*%zeta + x%*%beta
    	}
    } else if (!is.null(x)) {
    		Xb <- x%*%beta
    } else {
    		Xb <-0
    }
    z <- matrix(ncol = k - 1, nrow = length(y))
    for (i in 1:(k - 1)) {
    	z[, i] <- alpha[i] + Xb
    }
	grad<-matrix(0, nrow=length(y),ncol=length(theta))
	if (link=="logit") {
		for (j in 1:(k-1)){
			if (j==1){
				grad[,j]<-Ymat[,1]/(1+exp(z[,1]))- Ymat[,2]*(1+exp(z[,2]))*exp(z[,1])/((exp(z[,2])-exp(z[,1]))*(1+exp(z[,1])))
      		} else if (j < k-1 ){
				grad[,j]<-Ymat[,j]*(1+exp(z[,j-1]))*exp(z[,j])/((exp(z[,j])-exp(z[,j-1]))*(1+exp(z[,j])))-Ymat[,j+1]*(1+exp(z[,j+1]))*exp(z[,j])/((exp(z[,j+1])-exp(z[,j]))*(1+exp(z[,j])))
			} else if (j == k-1) {
				grad[,j]<- Ymat[,k-1]*(1+exp(z[,k-2]))*exp(z[,k-1])/((exp(z[,k-1])-exp(z[,k-2]))*(1+exp(z[,k-1])))-Ymat[,k]*exp(z[,k-1])/(1+exp(z[,k-1]))
			}
		}	
	} else if (link=="probit") {
		for (j in 1:(k-1)){
			if (j==1){
				grad[,j]<-Ymat[,1]*dnorm(z[,1])/pnorm(z[,1])- Ymat[,2]*(dnorm(z[,1])/(pnorm(z[,2])-pnorm(z[,1])))
			} 
			else if (j < k-1 ){
				grad[,j]<-Ymat[,j]*(dnorm(z[,j])/(pnorm(z[,j])-pnorm(z[,j-1])))-Ymat[,j+1]*(dnorm(z[,j])/(pnorm(z[,j+1])-pnorm(z[,j])))
			}
			else if (j == k-1) {
				grad[,j]<- Ymat[,k-1]*(dnorm(z[,k-1])/(pnorm(z[,k-1])-pnorm(z[,k-2])))-Ymat[,k]*dnorm(z[,k-1])/(1-pnorm(z[,k-1]))
			}
		}
	} else if (link=="cloglog") {
		for (j in 1:(k-1)){
			if (j==1){
				grad[,j]<-(t(Ymat*((exp(z[,j]))/(exp(exp(z[,j]))-1))))[j,]
			} else if (j < k-1 ){
				grad[,j]<-(t(Ymat*((-exp(exp(z[,j])))*(exp(exp(z[,j])+z[,j])-exp(z[,j])))))[j,]
			} else if (j == k-1) {
				grad[,j]<- (t(Ymat*(-exp(z[,j]))))[j,]
			}
		}
	}
	c(-apply(grad,2,sum))
}