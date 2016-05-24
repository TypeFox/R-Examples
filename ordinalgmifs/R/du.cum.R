du.cum <-
function(z, Ymat, k, x, link) {
	if (link=="logit") {
       	u <- matrix(ncol = k, nrow = dim(x)[2])
       	u.a <- u.a2 <- numeric()	
        for (j in 1:k) {
            # derivative for 1st logit
            if (j == 1) {
                mult <- 1 / (1 + exp(z[,1]))
            }
            # derivative for 2 to k-1 logits
            else if (j <= k-1) {
                mult <- -(exp(z[,j] + z[,j-1]) - 1) / 
                        ((1 + exp(z[,j])) * (1 + exp(z[,j-1])))
            }
            # derivative for kth logit 
            else if (j == k) {
                mult <- -exp(z[,j-1]) / (1 + exp(z[,j-1]))
            }
            u[,j] <- (t(Ymat * mult) %*% x)[j, ]
        }
        update.value <- min(rowSums(-u, na.rm=TRUE))
        update.j <- which.min(rowSums(-u, na.rm=TRUE))  
    } else if (link=="probit") {
    	u1<-matrix(nrow = dim(x)[1],ncol = k)
		for (j in 1:k){
			if (j==1){
				u1[,j]<-Ymat[,1]*dnorm(z[,1])/pnorm(z[,1])
			} else if (j <= k-1 ){
				u1[,j]<-Ymat[,j]*(dnorm(z[,j])-dnorm(z[,j-1]))/(pnorm(z[,j])-pnorm(z[,j-1]))
			} else if (j == k) {
				u1[,j]<- -Ymat[,k]*dnorm(z[,k-1])/(1-pnorm(z[,k-1])+1e-16)
			}
		}
		u<- -t(x) %*% rowSums(u1, na.rm=TRUE)
    	update.value <- min(u)
    	update.j <- which.min(u)
    } else if (link=="cloglog") {
		u <- matrix(ncol = k, nrow = dim(x)[2])
 		for (j in 1:k) {
		# derivative for 1st cloglog
			if (j==1) {
				mult<-(exp(z[,j]))/(exp(exp(z[,j]))-1)
			# derivative for 2 to k-1 cloglogs
			} else if (j <=k-1) {
				mult<- -(exp(exp(z[,j-1])+z[,j])-exp(exp(z[,j])+z[,j-1]))/(exp(exp(z[,j-1]))-exp(exp(z[,j])))
			# derivative for kth cloglog 
			} else if (j==k) {
				mult<- -exp(z[,j-1])
			}
			u[,j] <- (t(Ymat * mult) %*% x)[j, ]     
		}
		update.value <- min(rowSums(-u, na.rm=TRUE))
        update.j <- which.min(rowSums(-u, na.rm=TRUE))  
    }
    list(update.j=update.j, update.value=update.value)  
}
