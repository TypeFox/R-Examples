##########################  Function part ##################

library(MASS)

##Function for AUC##
AUC <- function(beta, Z, lambda){ 
		a<-sum(1/(1+exp(-lambda*t(beta)%*%t(Z))))/nrow(Z)
		a
}

## Function for AUC when input is X and Y
# Y: diseased sample
# X: nondiseased sample
AUC.emp <- function(X, Y){
  mean(outer(Y,X, FUN=function(y,x) 1*(y>x)))
}


##function for hessian matrix of AUC##
hessAUC <- function(beta, X, Y, lambda){
		m<-nrow(X)
		n<-nrow(Y)
		trait<-ncol(X)
		a<-matrix(0, trait, trait)
			for (i in 1 : m){
				for (j in 1 : n){
					exp.tmp<-exp(-lambda*t(beta)%*%(Y[j,]-X[i,]))
					a<-lambda^2*(Y[j,]-X[i,])%*%t(Y[j,]-X[i,])*c(exp.tmp*(exp.tmp-1)/(1+exp.tmp)^3)+a
					}
				}
		a<- a/(m*n)
		a
}

### The function of grad_square in the GCV ###

gradsqr <- function(beta, X, Y, lambda){
    m<-nrow(X)
    n<-nrow(Y)
    trait<-ncol(X)
    a<-a1<-a2<-matrix(0, trait, trait)
    
    for (i in 1:m){
      an<-rep(0, trait)
      for (j in 1:n){
        grad <- lambda*(Y[j,]-X[i,])*exp(-lambda*t(beta)%*%(Y[j,]-X[i,]))/(1+exp(-lambda*t(beta)%*%(Y[j,]-X[i,])))^2
        an <- grad + an
      }
      a1<-an%*%t(an)/n^2+a1
    }
    
    a1<-a1/m^2
    
    for (j in 1:n){
      am<-rep(0, trait)
      for (i in 1:m){
        grad <- lambda*(Y[j,]-X[i,])*exp(-lambda*t(beta)%*%(Y[j,]-X[i,]))/(1+exp(-lambda*t(beta)%*%(Y[j,]-X[i,])))^2
        am <- grad + am
      }
      a2<-am%*%t(am)/m^2+a2
    }
    
    a2<-a2/n^2
    
    for (i in 1 : m){
      for (j in 1 : n){
        grad<-lambda*(Y[j,]-X[i,])*exp(-lambda*t(beta)%*%(Y[j,]-X[i,]))/(1+exp(-lambda*t(beta)%*%(Y[j,]-X[i,])))^2
        a<-grad%*%t(grad)+a
      }
    }
    a <- a/(m*n)^2
    a0 <- a1+a2-a
    a0
}


### Function for gradient of AUC after applying Lagrange Multiplyer ###
gradAUC.Lang <- function(par, Z, lambda){ 
		a<-rep(0,ncol(Z))
		beta<-matrix(par[1:ncol(Z)],ncol(Z),1)
		a<-colSums(Z*lambda*c(exp(-lambda*t(beta)%*%t(Z))/(1+exp(-lambda*t(beta)%*%t(Z)))^2))/nrow(Z) + 2*par[ncol(Z)+1]*t(beta)
		a0 <- sum(beta^2)-1
		g <- c(a,a0)
		g
}

###  Function for estimating \beta using kernal function ###

betahat <- function(X, Y, init, lambda){
	m<-nrow(X)
	n<-nrow(Y)
	Z<-matrix(0, m*n, ncol(X))
	for(i in 1:m){
		for(j in 1:n){
			Z[(i-1)*n+j, ] <- (Y[j, ]-X[i, ])
		}
	}
	# estimate beta via Lagrange Multiplyer #
	result <- nlsolve(c(init,-0.1), gradAUC.Lang, Z=Z, lambda=lambda)
	beta.hat <- result$ans[1:ncol(Z)]
	converge <- result$converge

	list(beta.hat=beta.hat, converge=converge)
}


### Function to translate beta into theta, the n-sphere constrain ###
beta2theta <- function(beta){
	dim <- length(beta)
	theta<-rep(0,dim-1)
	for (i in 1:dim-1){
		theta[i]<-atan(sqrt(sum(beta[c((i+1):dim)]^2))/beta[i])
	}
	theta	
}

### Function to translate theta to beta, the n-sphere constrain ###
theta2beta <- function(theta){
		dim<-length(theta)+1
		beta<-matrix(1,dim,1)
### assign the n-sphere coordination #######
		for (k in 1: dim){
			if(k==1){
				beta[1,1]<-cos(theta[1]);				}else if (k== dim){
					for(ii in 1:(k-1)){
							beta[k,1]<-beta[k,1]*sin(theta[ii])					}
			}else {
					for(ii in 1:(k-1)){
						beta[k,1]<-beta[k,1]*sin(theta[ii])					}
					beta[k,1]<-beta[k,1]*cos(theta[k])			}		
		}		
	beta
}


##########################################
nlsolve <- function(par, fn, method="BFGS", nstarts = 1, ...) {
###########################
func <- function(x, ...) sum(fn(x, ...)^2)
ans <- optim(par, fn=func, method=method, ...)
err <- sqrt(ans$val/length(par)) 
istart <- 1
index<-1;

while (err > 0.01 & istart < nstarts) {
par <- par * ( 1 + runif(length(par), -1, 1) )  # generating random starting values
ans <- try (optim(par, fn=func, method=method, ...))
if (class(ans) != "try-error") err <- sqrt(ans$val/length(par)) 
istart <- istart + 1
}

#if (err > 0.01) cat(" \n Caution:  Zero of the function has not been located!!! \n Increase nstart \n \n")
# if not converge, try the global search
	if (err > 0.001) { 
		index<-0;
		
	}

list(ans=ans$par, converge=index)
}


optAUC <- function(X, Y, column.select=c(1:ncol(X)), lambda=5, scale=TRUE){
	cross<-0
	mdim<-nrow(X)
	ndim<-nrow(Y)
  if(scale==TRUE){
    all.data <- scale(rbind(X,Y))
	  cov <- var(all.data)
    X <- all.data[c(1:mdim),]
    Y <- all.data[-c(1:mdim),]
    Y.mean <- colSums(Y)/nrow(Y)
    X.mean <- colSums(X)/nrow(X)
  }else{
    all.data <- scale(rbind(X,Y))
    cov <- var(all.data)
    Y.mean <- colSums(Y)/nrow(Y)
    X.mean <- colSums(X)/nrow(X)
  }
	if (is.matrix(X[,column.select])){
		XC<-X[,column.select]
		YC<-Y[,column.select]
		int0<-solve(cov[column.select,column.select]+cov[column.select,column.select])%*%(Y.mean[column.select]-X.mean[column.select])
		init<-int0/sqrt(sum(int0^2))
		result<-betahat(XC,YC,init,lambda)
		beta<-result$beta.hat
		converge <- result$converge

		for(i in 1:mdim){
			for(j in 1:ndim){
				u0<-t(beta)%*%(YC[j,]-XC[i,]);
				cross<-1/(1+exp(-lambda*u0))+cross
			}
		}
		correct <- sum(diag(solve(hessAUC(beta, XC, YC, lambda))%*%gradsqr(beta, XC, YC, lambda)))
		
	} else {
		converge<-1
		int0<-solve(cov[column.select,column.select]+cov[column.select,column.select])%*%(Y.mean[column.select]-X.mean[column.select])
		init<-int0/sqrt(sum(int0^2))
		beta<-init
		XC<-t(t(X[,column.select]))
		YC<-t(t(Y[,column.select]))
		for(i in 1:mdim){
			for(j in 1:ndim){
				u0<-YC[j,]-XC[i,]
				cross<-1/(1+exp(-lambda*u0))+cross
			}
		}
		correct <- sum(diag(solve(hessAUC(beta, XC, YC, lambda))%*%gradsqr(beta, XC, YC, lambda)))     
	}
	CV <- cross/(mdim*ndim)
	GCV <- CV - abs(correct)
	list(ACV=CV, GCV=GCV, beta=beta, converge=converge)
}
