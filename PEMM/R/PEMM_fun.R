PEMM_fun <-
function(X, phi, lambda=NULL, K=NULL,   pos=NULL, tol=0.001, maxIter=100){


        
	get.llik <- function(X, mu, S, phi, phi0){
	    if (length(which(X<0))==0) pos<- TRUE else pos<- FALSE
	    
	    llik <- 0
	    for (j in 1:nrow(X)){
	      Xi <- X[j,]
	      idxi <- which(!is.na(Xi))
	      Xi <- Xi[idxi]
	      Si <- as.matrix(S[idxi,idxi])
	      Sii <- my.solve(Si)
	      oo <-    - log(det(Si))  -(Xi-mu[idxi])%*%Sii%*%(Xi-mu[idxi]) 
	      oo <- 0.5*oo
	      nmis <- length(X[j,])-length(idxi)
	      if (phi==0){
	          vv <- 0
	      } else {
	         if (pos){
	
		  vv1=sum(log(1-exp(-phi0-phi*Xi)))
		  vv2 <- -phi0*nmis
		  if (nmis==0){
		     vv3 <- 0
		  } else {
		      Smi <- as.matrix(S[-idxi,-idxi])
		      mu.mis <- mu[-idxi]+matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%(Xi-mu[idxi])
		      S.mis <- Smi - matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%matrix(S[idxi,-idxi],ncol=nmis)
		      vv3<- -sum(mu.mis)*phi+0.5*sum(S.mis)*phi^2  
		  }
	          vv <- vv1+vv2+vv3
	        } else {
	            vv1=sum(log(1-exp(-phi0-phi*Xi^2)))
		    vv2 <- -phi0*nmis
		    if (nmis==0){
		  	     vv3 <- 0
		    } else {
		  	  Smi <- as.matrix(S[-idxi,-idxi])
		  	  mu.mis <- mu[-idxi]+matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%(Xi-mu[idxi])
		  	  S.mis <- Smi - matrix(S[-idxi,idxi],nrow=nmis)%*%Sii%*%matrix(S[idxi,-idxi],ncol=nmis)
		  	  Smis.inv <- my.solve(S.mis)
		  	  Smii <- Smis.inv
		  	  diag(Smis.inv) <- diag(Smis.inv)+2*phi
		  	  A <- my.solve(Smis.inv)
	  	  
		  	  vv3<- 0.5*(log(det(A)) -log(det(S.mis)) + matrix(mu.mis,nrow=1)%*%(Smii%*%A%*%Smii-Smii)%*%matrix(mu.mis,ncol=1) ) 
		    }
	            vv <- vv1+vv2+vv3
	
	        }
	      }
	      llik <- llik+ oo+vv
	    }  
    
	    pllik <- llik
	    return(llik) 
	}

	my.solve <- function(X){
	   if (!is.matrix(X))  X <- matrix(X, nrow=sqrt(length(X)))
	   ss <- svd(X)
	   Xinv <- ss$u%*%diag(1/ss$d, nrow=nrow(X), ncol=nrow(X))%*%t(ss$v)
	   return(Xinv)
	}

	gets1 <- function(sigma, phi){
	  p <- nrow(sigma)
	  s1 <- my.solve(sigma)+diag(2*phi, p, p)
	  return(s1)
	}

	get.bg <- function(sigma, mu, phi){
	    p <- length(mu)
	    s1 <- gets1(sigma, phi)
	    s1inv <- my.solve(s1)
    
	    ccen <- my.solve(sigma)
	    A <- s1inv%*%ccen
	    beta <- A%*%mu

	    gamma <- s1inv 
    
	    return(list(beta=beta, gamma=gamma))
	}

     	find.lambda <- function(Sigma,N,p,K, delta=delta){
	    ffL <- function(lambda, Sigma, N, p, K){
	        Sigma.new <- N/(N+K)*Sigma*(N-1)/N + lambda/(N+K)*diag(1, p, p)
	        return(abs(min(as.double(eigen(N*Sigma.new)$value))))
	    }

	    Sigma2 = Sigma
	    while (!is.double(eigen(N*Sigma2)$value)){
	        delta = delta+1
	        Sigma2 = N/(N+K)*Sigma*(N-1)/N + delta/(N+K)*diag(1, p, p)
	    }
	    Sigma=Sigma2
	
	    oo <- -min(as.double(eigen(N*Sigma)$value))
	    if (oo>0){
	        lambda <- optimize(ffL, lower=0,upper=oo, Sigma=Sigma, N=N, p=p, K=K)$minimum+delta
	    } else {
	        lambda <- delta
	    }
	    return(lambda)
       }

       if (is.null(pos)){
         if (length(which(X<0))==0) pos<- TRUE else pos<- FALSE
       }

          if (phi==0) {
           phi0 <- -log(mean(is.na(X))) 
       } else {
           phi0 <- 0
       }
       
       p <- ncol(X)
       N <- nrow(X)
       if (is.null(K)){
         K=5
       } 
       X.hat <- X
       
       ## initial estimates
       mu.new <- matrix(colMeans(X,na.rm=T),ncol=1)
       Sigma.new <- cov(X, use="pairwise.complete")
       Sigma.new[is.na(Sigma.new)] <- 0  

       diff <- 999
       iter <- 0
       if (is.null(lambda)) {
         delta=5
         Lambda <- find.lambda(Sigma.new,N=N,p=p,K=K, delta=delta)
       } else {
         Lambda <- lambda
       }
       Sigma.new <- N/(N+K)*Sigma.new*(N-1)/N + Lambda/(N+K)*diag(1, p, p)
       illik <- 999
       while(iter<maxIter & diff>tol){
         iter <- iter+1
         mu <- mu.new
         Sigma <- Sigma.new
         
         cov.add <- matrix(0, p, p)
         for (i in 1:nrow(X)){
           ii <- which(is.na(X[i,]))
           if (length(ii)>=1){
             Soo <- as.matrix(Sigma[-ii,-ii])
             pi <- nrow(Soo)
             mu.mis <-  mu[ii]+Sigma[ii,-ii]%*%my.solve(Soo)%*%(X[i,-ii]-mu[-ii]) 
             mu.mis <- matrix(mu.mis,ncol=1)
             cov.mis <- Sigma[ii,ii] - Sigma[ii,-ii]%*%my.solve(Soo)%*%Sigma[-ii, ii]
             
             if (phi!=0 & pos==TRUE){

                 X.hat[i,ii]<- mu.mis - phi*cov.mis%*%matrix(1,nrow=length(mu.mis),ncol=1)
                 cov.ii <- cov.mis
                
             } else if (phi!=0 & pos==FALSE){

                 oo <- get.bg(cov.mis, mu.mis, phi)             
                 X.hat[i,ii] <- oo$beta
                 cov.ii <- oo$gamma

             } else if (phi==0) {

                 X.hat[i,ii] <- mu.mis
                 cov.ii <- cov.mis 

             }     
             cov.add[ii,ii] <- cov.add[ii,ii] + cov.ii
           }
         }          
         
         mu.new <- colMeans(X.hat)
         Sig <- cov(X.hat)*(N-1)/N+cov.add/N
         if (is.null(lambda)) {
	      Lambda <- find.lambda(Sig,N=N,p=p,K=K,delta=delta)
	 } else {
	      Lambda <- lambda
         }
         Sigma.new <- N/(N+K)*(Sig)  + Lambda/(N+K)*diag(1, p, p)

         illik <- c(illik, get.llik(X, mu.new, Sigma.new,phi=phi,phi0=phi0) - sum(diag(Lambda*my.solve(Sigma.new))) -K*log(det(Sigma.new)) ) 
         diff <- abs(illik[iter+1]-illik[iter])/abs(illik[iter])
         
             if (is.na(diff)) {
            diff <- 0  
            warning("The algorithm does not converge!")
         }
       }

       if (iter==maxIter) warning("The algorithm does not converge!")
      
       return(list(mu=mu, Sigma=Sigma, Xhat=X.hat))
}

