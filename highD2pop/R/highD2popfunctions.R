###############################################################################################
###############################################################################################
################### DATA GENERATOR  ###########################################################
###############################################################################################
###############################################################################################

build2popData <- function(	n,
									m,
									p,
									muX,
									muY, 
									dep, 
									commoncov, 
									VarScaleY, 
									ARMAparms,
									LRparm,
									S=1, 
									innov=function(n,...) rnorm(n,0,1),
									heteroscedastic=FALSE,
									het.diag)
{
	
# # test parms:								
# n = 45
# m = 60
# p = 300
# muX = rep(0,p)
# muY = rep(0,p)
# dep = 'LR'
# commoncov = FALSE
# VarScaleY = 5
# ARMAparms = list(coefs=list(ar=c(.4,.3),ma=c(.3,.5)))
# LRparm = .55
# S=3
# innov=function(n,...) rnorm(n,0,1)
# heteroscedastic=TRUE
# het.diag = diag(1 + rexp(p))
	    
    # create "errors"
    
   	ERRORS <- list()
    
    if(dep == 'IND'){
    	
	        for(s in 1:S)
	        {
	
    		    Zx <- matrix(innovate(n*p,rand.gen=innov), nrow = n)
	            Zy <- matrix(innovate(m*p,rand.gen=innov), nrow = m)   
	            
	            ERRORS[[length(ERRORS)+1]] <- list(Zx = Zx , Zy = Zy)
    		
    		}
    		
    		
    } else if(dep == 'ARMA'){
    	
    	Zx <- matrix(0,n,p)
        Zy <- matrix(0,m,p)
        
        for(s in 1:S)
	        {
	        	for(i in 1:max(n,m))
	            {
	                if (i <= n) Zx[i,] <- arima.sim( model = ARMAparms$coefs , n = p , rand.gen = innov)
	                if (i <= m) Zy[i,] <- arima.sim( model = ARMAparms$coefs , n = p , rand.gen = innov)    
	            }	
	            
	            ERRORS[[length(ERRORS)+1]] <- list(Zx = Zx , Zy = Zy)

	         }   
	            	
    } else if(dep == 'LR'){
    	
    	
	        H <- .5*(2 - LRparm)  
	        # Note that this part need only be done once per p,alpha combination. It is slow.
	        
	        R<-matrix(0,p,p)
	        
	        for(i in 1:p)
	        {
	            for(j in 1:p)
	            {
	                k <- abs(i - j)
	                R[i,j] <- .5*( (k + 1)^(2*H) + ((k - 1)^(2))^H - 2*k^(2*H))
	            }
	        }
	        
	        diag(R)<-1
	        
	        U <- chol(R)
	        
	        
	        for(s in 1:S)
	        {
				
				Z0x <- innovate(n*p,rand.gen=innov)
	            Zx <- matrix( Z0x, nrow = n) %*% U 
	            Z0y <- innovate(m*p,rand.gen=innov)
	            Zy <- matrix( Z0y , nrow = m) %*% U
	         
	         	ERRORS[[length(ERRORS)+1]] <- list(Zx = Zx , Zy = Zy)
	            
			}

    }
    
    # impose heteroscedasticity on innovations if heteroscedastic == TRUE
    
    if(heteroscedastic == TRUE)
    {
    	
    	for(s in 1:S)
    	{
    		
    			ERRORS[[s]]$Zx <- ERRORS[[s]]$Zx  %*% het.diag
    			ERRORS[[s]]$Zy <- ERRORS[[s]]$Zy  %*% het.diag
    		
    	}
    	
    }
    
    # ts.plot(ERRORS[[1]]$Zx[1,])
    
    # scale innovations differently if commoncov == FALSE
    
    if(commoncov == FALSE)
    {
    	
    	for(s in 1:S)
    	{
    		 
    			ERRORS[[s]]$Zy <- sqrt(VarScaleY)*ERRORS[[s]]$Zy  
    		
    	}
    	
    }	

        
    # now build matrices of means and add to them the error matrices to create data
	
    MUX <- matrix(muX,nrow=n,ncol=p,byrow=T)
    MUY <- matrix(muY,nrow=m,ncol=p,byrow=T)
    

    DATA <- list()
    
    for(s in 1:S)
    {
    
            X <- MUX + ERRORS[[s]]$Zx 
	       	Y <- MUY + ERRORS[[s]]$Zy
	       	
	       	DATA[[length(DATA)+1]] <- list( X = X,
	       									Y = Y,
	       									n = n,
											m = m,
											p = p,
											muX = muX,
											muY = muY, 
											dep = dep, 
											commoncov = commoncov, 
											VarScaleY = VarScaleY, 
											ARMAparms = ARMAparms,
											LRparm = LRparm,
											S = S, 
											innov=function(n,...) rnorm(n,0,1),
											heteroscedastic=heteroscedastic,
											het.diag = het.diag)
	       	
	} 
	

	
# ts.plot(ERRORS[[1]]$Zx[1,])
# lines(ERRORS[[1]]$Zy[1,],col="red")
       	
	return(DATA) 
        
}
    
    
build2popData_muYafunction <- function(	n,
										m,
										p,
										muX,
										muY, 
										dep, 
										commoncov, 
										VarScaleY, 
										ARMAparms,
										LRparm,
										S=1, 
										innov=function(n,...) rnorm(n,0,1),
										heteroscedastic=FALSE,
										het.diag)
{
	
# # test parms:								
# n = 45
# m = 60
# p = 300
# muX = rep(0,p)
# muY = rep(0,p)
# dep = 'LR'
# commoncov = FALSE
# VarScaleY = 5
# ARMAparms = list(coefs=list(ar=c(.4,.3),ma=c(.3,.5)))
# LRparm = .55
# S=3
# innov=function(n,...) rnorm(n,0,1)
# heteroscedastic=TRUE
# het.diag = diag(1 + rexp(p))
	    
    # create "errors"
    
   	ERRORS <- list()
    
    if(dep == 'IND'){
    	
	        for(s in 1:S)
	        {
	
    		    Zx <- matrix(innovate(n*p,rand.gen=innov), nrow = n)
	            Zy <- matrix(innovate(m*p,rand.gen=innov), nrow = m)   
	            
	            ERRORS[[length(ERRORS)+1]] <- list(Zx = Zx , Zy = Zy)
    		
    		}
    		
    		
    } else if(dep == 'ARMA'){
    	
    	Zx <- matrix(0,n,p)
        Zy <- matrix(0,m,p)
        
        for(s in 1:S)
	        {
	        	for(i in 1:max(n,m))
	            {
	                if (i <= n) Zx[i,] <- arima.sim( model = ARMAparms$coefs , n = p , rand.gen = innov)
	                if (i <= m) Zy[i,] <- arima.sim( model = ARMAparms$coefs , n = p , rand.gen = innov)    
	            }	
	            
	            ERRORS[[length(ERRORS)+1]] <- list(Zx = Zx , Zy = Zy)

	         }   
	            	
    } else if(dep == 'LR'){
    	
    	
	        H <- .5*(2 - LRparm)  
	        # Note that this part need only be done once per p,alpha combination. It is slow.
	        
	        R<-matrix(0,p,p)
	        
	        for(i in 1:p)
	        {
	            for(j in 1:p)
	            {
	                k <- abs(i - j)
	                R[i,j] <- .5*( (k + 1)^(2*H) + ((k - 1)^(2))^H - 2*k^(2*H))
	            }
	        }
	        
	        diag(R)<-1
	        
	        U <- chol(R)
	        
	        
	        for(s in 1:S)
	        {
				
				Z0x <- innovate(n*p,rand.gen=innov)
	            Zx <- matrix( Z0x, nrow = n) %*% U 
	            Z0y <- innovate(m*p,rand.gen=innov)
	            Zy <- matrix( Z0y , nrow = m) %*% U
	         
	         	ERRORS[[length(ERRORS)+1]] <- list(Zx = Zx , Zy = Zy)
	            
			}

    }
    
    # impose heteroscedasticity on innovations if heteroscedastic == TRUE
    
    if(heteroscedastic == TRUE)
    {
    	
    	for(s in 1:S)
    	{
    		
    			ERRORS[[s]]$Zx <- ERRORS[[s]]$Zx  %*% het.diag
    			ERRORS[[s]]$Zy <- ERRORS[[s]]$Zy  %*% het.diag
    		
    	}
    	
    }
    
    # ts.plot(ERRORS[[1]]$Zx[1,])
    
    # scale innovations differently if commoncov == FALSE
    
    if(commoncov == FALSE)
    {
    	
    	for(s in 1:S)
    	{
    		 
    			ERRORS[[s]]$Zy <- sqrt(VarScaleY)*ERRORS[[s]]$Zy  
    		
    	}
    	
    }	

        
    # now build matrices of means and add to them the error matrices to create data	
	
    MUX <- matrix(muX,nrow=n,ncol=p,byrow=T)
    
    MUY.list <- list()
    for(s in 1:S)
    {
    	MUY.list[[length(MUY.list)+1]] <- matrix(muY(),nrow=m,ncol=p,byrow=T)
    }

    DATA <- list()
    
    for(s in 1:S)
    {
    
            X <- MUX + ERRORS[[s]]$Zx 
	       	Y <- MUY.list[[s]] + ERRORS[[s]]$Zy
	       	
	       	DATA[[length(DATA)+1]] <- list( X = X,
	       									Y = Y,
	       									n = n,
											m = m,
											p = p,
											muX = muX,
											muY = MUY.list[[s]][1,], 
											dep = dep, 
											commoncov = commoncov, 
											VarScaleY = VarScaleY, 
											ARMAparms = ARMAparms,
											LRparm = LRparm,
											S = S, 
											innov=function(n,...) rnorm(n,0,1),
											heteroscedastic=heteroscedastic,
											het.diag = het.diag)
	       	
	} 
	

	
# ts.plot(ERRORS[[1]]$Zx[1,])
# lines(ERRORS[[1]]$Zy[1,],col="red")
       	
	return(DATA) 
        
}
    


###############################################################################################
###############################################################################################
################### INNOVATION-GENERATING FUNCTIONS ###########################################
###############################################################################################
###############################################################################################

innovate<-function (n, rand.gen = rnorm, innov = rand.gen(n, ...), ...) 
{
	    
    return(innov)
}



rgammashift <- function(n,shape,scale)
{
	x<- rgamma(n=n,shape=shape,scale=scale) - shape*scale
	return(x)
}


rdblepareto <- function(n,shape,scale)
{
# Note that the first moment is infinite if the shape parameter is less than or equal to 1.
	
	u<-runif(n,0,1)
	x <- scale*(1-u)^(-1/shape) - scale
	y <- sample(c(-1,1),n,replace=TRUE)*x
	return(y)
	
}


###############################################################################################
###############################################################################################
################### GCT TEST ##################################################################
###############################################################################################
###############################################################################################


center <- function(xy,n,m,ntoorderminus=2)
{
	
	if(ntoorderminus == 0) { return(1)
	} else if (ntoorderminus==1) {
	x <- xy[1:n]
	y <- xy[(n+1):(n+m)]
	
	sig.sq.x.hat <- mean(x^2) - mean(x)^2
	sig.sq.y.hat <- mean(y^2) - mean(y)^2
	
	tau.sq.hat <- sig.sq.x.hat + (n/m) * sig.sq.y.hat

	mu.3.hat  <- mean((x - mean(x))^3)
	eta.3.hat <- mean((y - mean(y))^3)

	mu.4.hat <-  mean((x - mean(x))^4)
	eta.4.hat <- mean((y - mean(y))^4)
	
	mu.5.hat <-  mean((x - mean(x))^5)
	eta.5.hat <- mean((y - mean(y))^5)
	
	
	a1 <- tau.sq.hat^(-1) * ( sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat )
	
	a2 <- tau.sq.hat^(-3) * 2 * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2
	

	c <- 1 +  n^(-1) * ( a1 + a2 ) 	

	return(c)	
	
	}	else if(ntoorderminus==2) {
	
	x <- xy[1:n]
	y <- xy[(n+1):(n+m)]
	
	sig.sq.x.hat <- (mean(x^2) - mean(x)^2)
	sig.sq.y.hat <- mean(y^2) - mean(y)^2
	
	tau.sq.hat <- sig.sq.x.hat + (n/m) * sig.sq.y.hat

	mu.3.hat  <- mean((x - mean(x))^3)
	eta.3.hat <- mean((y - mean(y))^3)

	mu.4.hat <-  mean((x - mean(x))^4)
	eta.4.hat <- mean((y - mean(y))^4)
	
	mu.5.hat <-  mean((x - mean(x))^5)
	eta.5.hat <- mean((y - mean(y))^5)
	
	
	a1 <- tau.sq.hat^(-1) * ( sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat )
	
	a2 <- tau.sq.hat^(-3) * 2 * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2
	
	
	b1 <- tau.sq.hat^(-2) * ( (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat ) - ((mu.4.hat - 3 * sig.sq.x.hat^2 ) + (n/m)^4 * (eta.4.hat - 3 * sig.sq.y.hat^2 )))
	
	b2 <- tau.sq.hat^(-3) * ( (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat ) * ( (mu.4.hat -  sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 ) ) - 4* ( mu.3.hat + (n/m)^2 * eta.3.hat) * ( mu.3.hat + (n/m)^3 * eta.3.hat) - 2 * ( mu.3.hat^2 + (n/m)^5 * eta.3.hat^2)  )
	
	b3 <- tau.sq.hat^(-4) * ( 6 * (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat ) * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2 - 6 * (mu.3.hat + (n/m)^2 * eta.3.hat) * ( mu.5.hat - 2 * mu.3.hat * sig.sq.x.hat +  (n/m)^4 * ( eta.5.hat - 2 * eta.3.hat * sig.sq.y.hat) ) - 3 * ((mu.4.hat - sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 ))^2 )
	
	b4 <- tau.sq.hat^(-5) *( 3 * (sig.sq.x.hat + (n/m) * sig.sq.y.hat ) * ((mu.4.hat - sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 ))^2 + 12 * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2  * ((mu.4.hat - sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 )) ) 
	
	c <- 1 +  n^(-1) * ( a1 + a2 ) + n^(-2) * ( b1 + b2 + b3 + b4)



	return(c)
	
	}

}



pzen <- function(r)
{
	x1 <- (1:floor((r-1)/2))/r
	y <- 1 - 6*x1^2 + 6*abs(x1)^3
	
	x2 <-  ((floor((r-1)/2)+1):r)/r	
	z <- 2*(1 - abs(x2))^3	
	
	vec <- c(y,z)	
    
	return(vec)		
}



trapezoid <- function(r)
{
    
	y <- rep(1,ceiling(r/2))
	
	s <-  ((ceiling(r/2)+1):r)	
	z <- 1 - (s - ceiling(r/2))/(r-ceiling(r/2))
	
	vec <- c(y,z)	
    
	return(vec)		
	
}

GCT.test <- function(X,Y,r,smoother="parzen",ntoorderminus=2)
{
	
	
	X.bar <- apply(X,2,mean)	
	Y.bar <- apply(Y,2,mean)	
	X.var <- apply(X,2,var)
	Y.var <- apply(Y,2,var)
	
	n <- length(X[,1])
	m <- length(Y[,1])
	p <- length(X[1,])
	
	D.sq <- (X.bar - Y.bar)^2 / (X.var / n  +  Y.var / m )
	
	T <- mean(D.sq)
	
	gamma <- acf(D.sq,lag.max=r,plot=FALSE,type="covariance")$acf
	
	if(smoother=="parzen"){	var.est <- (2* t(pzen(r)) %*% gamma[-1] + gamma[1]) / p
	} else if(smoother== "trapezoid") {var.est <- (2* t(trapezoid(r)) %*% gamma[-1] + gamma[1]) / p
        
        if(var.est < 0) 
        { 
            var.est <- abs(var.est)
            smoother <- "trapezoid (neg)"
        }
    }
    
    
	center.est <- apply(rbind(X,Y),2,center,n=n,m=m,ntoorderminus=ntoorderminus)
	
    TSvalue <- (T - mean(center.est)) / sqrt(var.est)
    
    pvalue <- 2*(1 - pnorm(abs(TSvalue),0,1))
    
    list(TSvalue=TSvalue, pvalue = pvalue, smoother= smoother, center=mean(center.est), T = T, var = var.est)
    
}


GCT.sim <- function(DATA,r,smoother="parzen",ntoorderminus=2)
{
	
	
	S <- length(DATA)
	pvalues <- TSvalues <- TT <- Var <- Cent <- numeric(S)
	smoothtype <- rep(smoother,S)			
    
	n <- DATA[[1]]$n
	m <- DATA[[1]]$m
	p <- DATA[[1]]$p
    
	if(smoother=="parzen")
	{
		for(s in 1:S)
		{
			X <- DATA[[s]]$X
			Y <- DATA[[s]]$Y
            
			X.bar <- apply(X,2,mean)	
			Y.bar <- apply(Y,2,mean)	
			X.var <- apply(X,2,var)
			Y.var <- apply(Y,2,var)
            
			D.sq <- (X.bar - Y.bar)^2 / (X.var / n  +  Y.var / m )
            
			T <- mean(D.sq)
            
			gamma <- acf(D.sq,lag.max=r,plot=FALSE,type="covariance")$acf
            
			var.est <- (2 * t(pzen(r)) %*% gamma[-1] + gamma[1]) / p 
            
			center.est <- apply(rbind(X,Y),2,center,n=n,m=m,ntoorderminus=ntoorderminus)
		    
			TSvalue <- (T - mean(center.est)) / sqrt(var.est)
            
		 	TSvalues[s] <- TSvalue
		 	
		 	pvalues[s] <-2*(1 - pnorm(abs(TSvalue),0,1))
		 	
		 	TT[s] <- T
		 	Cent[s] <- mean(center.est)
		 	Var[s] <- var.est
		 	
        }
	}
    
	if(smoother=="trapezoid")
	{
		for(s in 1:S)
		{
			X <- DATA[[s]]$X
			Y <- DATA[[s]]$Y
            
			X.bar <- apply(X,2,mean)	
			Y.bar <- apply(Y,2,mean)	
			X.var <- apply(X,2,var)
			Y.var <- apply(Y,2,var)
            
			D.sq <- (X.bar - Y.bar)^2 / (X.var / n  +  Y.var / m )
            
			T <- mean(D.sq)
            
			gamma <- acf(D.sq,lag.max=r,plot=FALSE,type="covariance")$acf
			
			var.est <- (2* t(trapezoid(r)) %*% gamma[-1] + gamma[1]) / p
            
	 		if(var.est < 0) 
  			{ 
				var.est <- abs(var.est)
	 			smoothtype[s] <- "trapezoid (neg)"
 			}
            
			center.est <- apply(rbind(X,Y),2,center,n=n,m=m,ntoorderminus=ntoorderminus)
		    
			TSvalue <- (T - mean(center.est)) / sqrt(var.est)
            
		 	TSvalues[s] <- TSvalue
		 	
		 	pvalues[s] <-2*(1 - pnorm(abs(TSvalue),0,1))
		 	
		 	TT[s] <- T
		 	Cent[s] <- mean(center.est)
		 	Var[s] <- var.est
		 	
        }
	}
    
    
	results <- data.frame( TSvalues=TSvalues , pvalues = pvalues , smoothtype=smoothtype , T = TT , cent = Cent, var = Var)
    
	return(results)
    
}


nonmiss <- function(x)
{
	
	return(-sum(is.na(x)-1))

}

center.missing <- function(xy,n.orig,m.orig,ntoorderminus=2)
{
	
	if(ntoorderminus == 0) { return(1)
	} else if (ntoorderminus==1) {
		
	x.orig <- xy[1:n.orig]
	y.orig <- xy[(n.orig+1):(n.orig+m.orig)]
	
	n <- sum(is.na(x.orig)==FALSE)
	m <- sum(is.na(y.orig)==FALSE)
	
	x <- x.orig[is.na(x.orig)==FALSE]
	y <- y.orig[is.na(y.orig)==FALSE]

	sig.sq.x.hat <- mean(x^2) - mean(x)^2
	sig.sq.y.hat <- mean(y^2) - mean(y)^2
	
	tau.sq.hat <- sig.sq.x.hat + (n/m) * sig.sq.y.hat

	mu.3.hat  <- mean((x - mean(x))^3)
	eta.3.hat <- mean((y - mean(y))^3)

	mu.4.hat <-  mean((x - mean(x))^4)
	eta.4.hat <- mean((y - mean(y))^4)
	
	mu.5.hat <-  mean((x - mean(x))^5)
	eta.5.hat <- mean((y - mean(y))^5)
	
	a1 <- tau.sq.hat^(-1) * ( sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat )
	
	a2 <- tau.sq.hat^(-3) * 2 * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2
	
	c <- 1 +  n^(-1) * ( a1 + a2 ) 	

	return(c)	
	
	}	else if(ntoorderminus==2) {
	
	x.orig <- xy[1:n.orig]
	y.orig <- xy[(n.orig+1):(n.orig+m.orig)]
	
	n <- sum(is.na(x.orig)==FALSE)
	m <- sum(is.na(y.orig)==FALSE)
	
	x <- x.orig[is.na(x.orig)==FALSE]
	y <- y.orig[is.na(y.orig)==FALSE]
	
	sig.sq.x.hat <- mean(x^2) - mean(x)^2
	sig.sq.y.hat <- mean(y^2) - mean(y)^2
	
	tau.sq.hat <- sig.sq.x.hat + (n/m) * sig.sq.y.hat

	mu.3.hat  <- mean((x - mean(x))^3)
	eta.3.hat <- mean((y - mean(y))^3)

	mu.4.hat <-  mean((x - mean(x))^4)
	eta.4.hat <- mean((y - mean(y))^4)
	
	mu.5.hat <-  mean((x - mean(x))^5)
	eta.5.hat <- mean((y - mean(y))^5)
	
	
	a1 <- tau.sq.hat^(-1) * ( sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat )
	
	a2 <- tau.sq.hat^(-3) * 2 * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2
	
	
	b1 <- tau.sq.hat^(-2) * ( (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat ) - ((mu.4.hat - 3 * sig.sq.x.hat^2 ) + (n/m)^4 * (eta.4.hat - 3 * sig.sq.y.hat^2 )))
	
	b2 <- tau.sq.hat^(-3) * ( (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat ) * ( (mu.4.hat -  sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 ) ) - 4* ( mu.3.hat + (n/m)^2 * eta.3.hat) * ( mu.3.hat + (n/m)^3 * eta.3.hat) - 2 * ( mu.3.hat^2 + (n/m)^5 * eta.3.hat^2)  )
	
	b3 <- tau.sq.hat^(-4) * ( 6 * (sig.sq.x.hat + (n/m)^2 * sig.sq.y.hat ) * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2 - 6 * (mu.3.hat + (n/m)^2 * eta.3.hat) * ( mu.5.hat - 2 * mu.3.hat * sig.sq.x.hat +  (n/m)^4 * ( eta.5.hat - 2 * eta.3.hat * sig.sq.y.hat) ) - 3 * ((mu.4.hat - sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 ))^2 )
	
	b4 <- tau.sq.hat^(-5) *( 3 * (sig.sq.x.hat + (n/m) * sig.sq.y.hat ) * ((mu.4.hat - sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 ))^2 + 12 * ( mu.3.hat  + (n/m)^2 *eta.3.hat )^2  * ((mu.4.hat - sig.sq.x.hat^2 ) + (n/m)^3 * (eta.4.hat - sig.sq.y.hat^2 )) ) 
	
	c <- 1 +  n^(-1) * ( a1 + a2 ) + n^(-2) * ( b1 + b2 + b3 + b4)


	return(c)
	
	}

}


GCT.test.missing <- function(X,Y,r,smoother="parzen",ntoorderminus=2)
{
	

	n <- apply(X,2,nonmiss)
	m <- apply(Y,2,nonmiss)
	p <- length(X[1,])
	
	overallpctmiss <- sum((dim(X)[1] - n + dim(Y)[1] - m )) / (sum(m)+sum(n))
	
	pctmissperX <- (dim(X)[1] - n)/dim(X)[1]
	pctmissperY <- (dim(Y)[1] - m)/dim(Y)[1]
	
	remove <- ( n < 3 ) | ( m < 3)
	
	if(sum(remove)>0)
		{
		
				X <- X[,as.logical(1-remove)]
				Y <- Y[,as.logical(1-remove)]
				n <- n[as.logical(1-remove)]
				m <- m[as.logical(1-remove)]
				p <- p - sum(remove)
		}
	
	X.bar <- apply(X,2,mean,na.rm=TRUE)	
	Y.bar <- apply(Y,2,mean,na.rm=TRUE)	
	X.var <- apply(X,2,var,na.rm=TRUE)
	Y.var <- apply(Y,2,var,na.rm=TRUE)
			
	D.sq <- (X.bar - Y.bar)^2 / (X.var / n  +  Y.var / m )
	
	T <- mean(D.sq)
	
	gamma <- acf(D.sq,lag.max=r,plot=FALSE,type="covariance")$acf
	
	if(smoother=="parzen"){	var.est <- (2* t(pzen(r)) %*% gamma[-1] + gamma[1]) / p
	} else if(smoother== "trapezoid") {var.est <- (2* t(trapezoid(r)) %*% gamma[-1] + gamma[1]) / p
        
        if(var.est < 0) 
        { 
            var.est <- abs(var.est)
            smoother <- "trapezoid (neg)"
        }
    }
    
	center.est <- apply(rbind(X,Y),2,center.missing,n.orig=dim(X)[1],m.orig=dim(Y)[1],ntoorderminus=ntoorderminus)
	
    TSvalue <- (T - mean(center.est)) / sqrt(var.est)
    
    pvalue <- 2*(1 - pnorm(abs(TSvalue),0,1))
    
    list(TSvalue=TSvalue, pvalue = pvalue, smoother= smoother, center=mean(center.est), T = T, var = var.est, p = p , r = r , 
    		overallpctmiss = overallpctmiss, pctmissperX = pctmissperX, pctmissperY = pctmissperY )
    
}


###############################################################################################
###############################################################################################
################### CAI, LIU, AND XIA 2011 TEST  FOR EQUAL COVARIANCE MATRICES ################
###############################################################################################
###############################################################################################

CLX.Covtest<- function(X,Y)
{
	
	n1 <- dim(X)[1]
	n2 <- dim(Y)[1]
	p <- dim(X)[2]
	
	W1 <- X - matrix(1/n1,n1,n1) %*% X
	W2 <- Y - matrix(1/n2,n2,n2) %*% Y
	
	S1 <- t(W1) %*% W1 / n1
	S2 <- t(W2) %*% W2 / n2
	
	Theta1 <- Theta2 <- matrix(0,p,p)
	
	for(i in 1:n1)
	{
		Theta1 <- Theta1 + (1/n1) * (  W1[i,] %*% t(W1[i,]) - S1 )^2
	}
	
	for(i in 1:n2)
	{
		Theta2 <- Theta2 + (1/n1) * (  W2[i,] %*% t(W2[i,]) - S2 )^2
	}
	
	W <- ( S1 - S2 ) / sqrt(Theta1/n1 + Theta2/n2)
	
	M <- W^2
	
	M.n <- max(M)
	# qqnorm(W)

	TSvalue <- M.n - 4*log(p) + log(log(p))
	
	pvalue <- 1 - exp( - 1/sqrt(8*pi) * exp(-TSvalue/2))
	
	output <- list( TSvalue = TSvalue,
					pvalue = pvalue)
					
	return(output)
	
}


###############################################################################################
###############################################################################################
################### CAI, LIU, AND XIA 2013 TEST  FOR EQUAL MEAN VECTORS #######################
###############################################################################################
###############################################################################################

CLX.test.equalcov <- function(X,Y)
{
	
	X.bar <- apply(X,2,mean)
	Y.bar <- apply(Y,2,mean)
	n1 <- nrow(X)
	n2 <- nrow(Y)
	p <- ncol(X)
	n <- n1 + n2 - 2
	
	S <- ((n1-1)*cov(X) + (n2-1)*cov(Y)) / (n1 + n2)
	
	fast.out <- fastclime(S)
	
	Omega.hat <- fast.out$icov[[length(fast.out$icov)]]
	
	Z <- Omega.hat %*% (X.bar - Y.bar )
	
	X.omega.hat <- X %*% Omega.hat
	Y.omega.hat <- Y %*% Omega.hat
	
	Omega.hat.0 <- ((n1-1)*cov(X.omega.hat) + (n2-1)*cov(Y.omega.hat)) / ( n1 + n2 ) 
	
	TSvalue <- (n1*n2)/(n1 + n2) * max(Z^2 / diag(Omega.hat.0) ) - 2 * log(p) + log( log(p)) 

	pvalue <- 1 - exp(-(1/sqrt(pi))*exp(-TSvalue/2))

	return(list( TSvalue = TSvalue , pvalue = pvalue))

}


CLX.test.unequalcov <- function(X,Y)
{
	
	X.bar <- apply(X,2,mean)
	Y.bar <- apply(Y,2,mean)
	n1 <- nrow(X)
	n2 <- nrow(Y)
	p <- ncol(X)
	n <- n1 + n2 - 2
	
	S1 <- (n1-1)*cov(X)/n1
	S2 <- (n2-1)*cov(Y)/n2
	
	S.unpooled <- ( S1 + (n1/n2)*S2 ) 
	fast.out <- fastclime(S.unpooled)
	
	Omega.hat <- fast.out$icov[[length(fast.out$icov)]]
	
	Z <- Omega.hat %*% (X.bar - Y.bar )
	
	X.omega.hat <- X %*% Omega.hat
	Y.omega.hat <- Y %*% Omega.hat
	
	Omega.hat.0 <- (n1-1)*cov(X.omega.hat)/n1 + (n2-1)*cov(Y.omega.hat)/n2 * (n1/n2)
	
	TSvalue <- n1 * max(Z^2 / diag(Omega.hat.0) ) - 2 * log(p) + log( log(p) ) 

	pvalue <- 1 - exp(-(1/sqrt(pi))*exp(-TSvalue/2))
	
	return(list( TSvalue = TSvalue , pvalue = pvalue))
}




CLX.sim.equalcov <- function(DATA)
{
	
	S <- length(DATA)
	pvalues <- TSvalues <- numeric(S)

	n1 <- DATA[[1]]$n
	n2 <- DATA[[1]]$m
	p <- DATA[[1]]$p
	n <- n1 + n2 - 2
    

		for(s in 1:S)
		{
			
			X <- DATA[[s]]$X
			Y <- DATA[[s]]$Y
			
			X.bar <- apply(X,2,mean)
			Y.bar <- apply(Y,2,mean)
			n1 <- nrow(X)
			n2 <- nrow(Y)
			n <- n1 + n2 - 2
			
			S <- ((n1-1)*cov(X) + (n2-1)*cov(Y)) / (n1 + n2)
			fast.out <- fastclime(S)
			
			Omega.hat <- fast.out$icov[[length(fast.out$icov)]]
			
			Z <- Omega.hat %*% (X.bar - Y.bar )
			
			X.omega.hat <- X %*% Omega.hat
			Y.omega.hat <- Y %*% Omega.hat
			
			Omega.hat.0 <- ((n1-1)*cov(X.omega.hat) + (n2-1)*cov(Y.omega.hat)) / ( n1 + n2 ) 
			
			TSvalue <- (n1*n2)/(n1 + n2) * max(Z^2 / diag(Omega.hat.0) ) - 2 * log(p) + log( log(p)) 
		
			pvalue <- 1 - exp(-(1/sqrt(pi))*exp(-TSvalue/2))
			
		 	TSvalues[s] <- TSvalue
		 	
		 	pvalues[s] <- pvalue
		 	
        }

	results <- data.frame( TSvalues = TSvalues , pvalues = pvalues )
    
	return(results)
    
}


CLX.sim.unequalcov <- function(DATA)
{
	
	S <- length(DATA)
	pvalues <- TSvalues <- numeric(S)

	n1 <- DATA[[1]]$n
	n2 <- DATA[[1]]$m
	p <- DATA[[1]]$p
	n <- n1 + n2 - 2
    

		for(s in 1:S)
		{
			
			X <- DATA[[s]]$X
			Y <- DATA[[s]]$Y
			
			X.bar <- apply(X,2,mean)
			Y.bar <- apply(Y,2,mean)
			n1 <- nrow(X)
			n2 <- nrow(Y)
			n <- n1 + n2 - 2
			
			S1 <- (n1-1)*cov(X)/n1
			S2 <- (n2-1)*cov(Y)/n2
		
			S.unpooled <- ( S1 + (n1/n2)*S2 ) 
			fast.out <- fastclime(S.unpooled)
			
			Omega.hat <- fast.out$icov[[length(fast.out$icov)]]
			
			Z <- Omega.hat %*% (X.bar - Y.bar )
			
			X.omega.hat <- X %*% Omega.hat
			Y.omega.hat <- Y %*% Omega.hat
			
			Omega.hat.0 <- (n1-1)*cov(X.omega.hat)/n1 + (n2-1)*cov(Y.omega.hat)/n2 * (n1/n2)
			
			TSvalue <- n1 * max(Z^2 / diag(Omega.hat.0) ) - 2 * log(p) + log( log(p) ) 
	
			pvalue <- 1 - exp(-(1/sqrt(pi))*exp(-TSvalue/2))
			
		 	TSvalues[s] <- TSvalue
		 	
		 	pvalues[s] <- pvalue
		 	
        }

	results <- data.frame( TSvalues = TSvalues , pvalues = pvalues )
    
	return(results)
    
}




CLX.sim.Covtest <- function(DATA)
{
	
	
	S <- length(DATA)
	pvalues <- TSvalues <- rejectequalcov <- numeric(S)

	n1 <- DATA[[1]]$n
	n2 <- DATA[[1]]$m
	p <- DATA[[1]]$p
	n <- n1 + n2 - 2
    

		for(s in 1:S)
		{
			
			X <- DATA[[s]]$X
			Y <- DATA[[s]]$Y
			
			X.bar <- apply(X,2,mean)
			Y.bar <- apply(Y,2,mean)
			n1 <- nrow(X)
			n2 <- nrow(Y)
			n <- n1 + n2 - 2
			
			S1 <- (n1-1)*cov(X)/n1
			S2 <- (n2-1)*cov(Y)/n2
	
			Covtest <- CLX.Covtest(X,Y)
			
			if( Covtest$pvalue < .05 )
			{
				S.unpooled <- ( S1 + (n1/n2)*S2 ) 
				fast.out <- fastclime(S.unpooled)
				rejectequalcov[s] <- 1
				
				Omega.hat <- fast.out$icov[[length(fast.out$icov)]]
				
				Z <- Omega.hat %*% (X.bar - Y.bar )
				
				X.omega.hat <- X %*% Omega.hat
				Y.omega.hat <- Y %*% Omega.hat
				
				Omega.hat.0 <- (n1-1)*cov(X.omega.hat)/n1 + (n2-1)*cov(Y.omega.hat)/n2 * (n1/n2)
				
				TSvalue <- n1 * max(Z^2 / diag(Omega.hat.0) ) - 2 * log(p) + log( log(p) ) 
		
				pvalue <- 1 - exp(-(1/sqrt(pi))*exp(-TSvalue/2))
					
			} else {
				
				S.pooled <- ( n1 * S1 + n2 * S2 ) / n
				fast.out <- fastclime(S.pooled)
				rejectequalcov[s] <- 0
			
				Omega.hat <- fast.out$icov[[length(fast.out$icov)]]
				
				Z <- Omega.hat %*% (X.bar - Y.bar )
				
				X.omega.hat <- X %*% Omega.hat
				Y.omega.hat <- Y %*% Omega.hat
				
				Omega.hat.0 <- ((n1-1)*cov(X.omega.hat) + (n2-1)*cov(Y.omega.hat)) / ( n1 + n2 ) 
				
				TSvalue <- (n1*n2)/(n1 + n2) * max(Z^2 / diag(Omega.hat.0) ) - 2 * log(p) + log( log(p) ) 
		
				pvalue <- 1 - exp(-(1/sqrt(pi))*exp(-TSvalue/2))
			
			}
			
		 	TSvalues[s] <- TSvalue
		 	
		 	pvalues[s] <- pvalue
		 	
        }

	results <- data.frame( TSvalues = TSvalues , pvalues = pvalues , rejectequalcov = rejectequalcov )
    
	return(results)
    
}


###############################################################################################
###############################################################################################
################### SRIVASTAVA AND KUBOKAWA 2013 TEST #########################################
###############################################################################################
###############################################################################################

SK.test <- function(X,Y)
{
		
	n1 <- nrow(X)
	n2 <- nrow(Y)
	p <- ncol(X)
	n <- n1 + n2 - 2
	
	X.bar <- apply(X,2,mean)
	Y.bar <- apply(Y,2,mean)
	
	B <- (X.bar %*% t(X.bar) -  X.bar %*% t(Y.bar) - Y.bar %*% t(X.bar) + Y.bar %*% t(Y.bar)) *(n1*n2/(n1+n2))
	
	S <- ((n1-1)*cov(X) + (n2-1)*cov(Y))/n
	
	D <- diag(diag(S))
	
	R <- solve(sqrt(D)) %*% S %*%  solve(sqrt(D))
	
	c.p.n <- 1 + sum(diag(R %*% R)/p^(3/2))
	
	TSvalue <- (sum(diag(B %*% solve(D))) - n * p / ( n - 2)) / sqrt( 2 * c.p.n * (sum(diag(R %*% R)) - n^(-1)*p^2) )
    
    pvalue <- 2*(1 - pnorm(abs(TSvalue),0,1))
    
    list(TSvalue=TSvalue, pvalue = pvalue)
    
}



SK.sim <- function(DATA)
{
	
	
	S <- length(DATA)
	pvalues <- TSvalues <- numeric(S)

	n1 <- DATA[[1]]$n
	n2 <- DATA[[1]]$m
	p <- DATA[[1]]$p
	n <- n1 + n2 - 2
    

		for(s in 1:S)
		{
			X <- DATA[[s]]$X
			Y <- DATA[[s]]$Y

			X.bar <- apply(X,2,mean)
			Y.bar <- apply(Y,2,mean)
			
			B <- (X.bar %*% t(X.bar) - 2 * X.bar %*% t(Y.bar) + Y.bar %*% t(Y.bar)) *(n1*n2/(n1+n2))
			
			S <- ((n1-1)*cov(X) + (n2-1)*cov(Y))/n
			
			D <- diag(diag(S))
			
			R <- solve(sqrt(D)) %*% S %*%  solve(sqrt(D))
			
			c.p.n <- 1 + sum(diag(R %*% R)/p^(3/2))
			
			TSvalue <- (sum(diag(B %*% solve(D))) - n * p / ( n - 2)) / sqrt( 2 * c.p.n * (sum(diag(R %*% R)) - n^(-1)*p^2) )
		    
		    pvalue <- 2*(1 - pnorm(abs(TSvalue),0,1))

            
		 	TSvalues[s] <- TSvalue
		 	
		 	pvalues[s] <-2*(1 - pnorm(abs(TSvalue),0,1))
		 	
        }

	results <- data.frame( TSvalues=TSvalues , pvalues = pvalues )
    
	return(results)
    
}


###############################################################################################
###############################################################################################
################### SRIVASTAVA TEST ###########################################################
###############################################################################################
###############################################################################################

# # 
# Sriv.test<-function(X1,X2)
# {

	# n1 <- length(X1[,1])
	# n2 <- length(X2[,1])
	# p <- length(X1[1,])

	# X1.bar <- ( 1 / n1 ) * t(X1) %*% matrix(1,n1,1)

	# H1 <- diag(rep(1,n1)) - (1 / n1) * matrix(1,n1,n1)

	# S1 <- (1/n1) * t(X1) %*% H1 %*% X1

	# X2.bar <- ( 1 / n2 ) * t(X2) %*% matrix(1,n2,1)

	# H2 <- diag(rep(1,n2)) - (1 / n2) * matrix(1,n2,n2)

	# S2 <- (1/n2) * t(X2) %*% H2 %*% X2

	# S <- ( n1*S1 + n2*S2 ) / (n1 + n2 - 2)

	# n <- n1 + n2 - 2

	# S.plus <- ginv(S)

	# D.plus <- t(X1.bar - X2.bar) %*% S.plus %*% (X1.bar - X2.bar)

	# F.plus <- (max(n,p) - min(n,p) + 1) / (n * min(p,n) ) * 1/( 1/n1 + 1/n2 ) * D.plus

	# b.hat <- ( n - 1 ) * ( n + 2) / n^2 * (sum(diag(S))/p)^2 / ( (1/p) * (sum(diag(S %*% S)) - (1/n)*sum(diag(S))^2 ))

	# c.p.n <- sqrt(( (p - n + 1) / ( p + 1) ))

	# Z.plus <- c.p.n * sqrt(n/2) * ( b.hat * F.plus - 1)

	# p.value_z <- 2*( 1 - pnorm(abs(Z.plus)) )

	# p.value_f <- 1 - pf(F.plus,min(n,p),max(n,p) - min(n,p) + 1)

	# return(list(F.plus = F.plus , p.value_f = p.value_f, Z.plus= Z.plus, p.value_z = p.value_z))

# }


# Sriv.sim<-function(DATA)
# {
	
	# S <- length(DATA)
	# p.values_f <- F.plus.values <- p.values_z <- Z.plus.values <- numeric(S)
				
	# n1 <- DATA[[1]]$n
	# n2 <- DATA[[1]]$m
	# p <- DATA[[1]]$p

	# for(s in 1:S)
	# {
		# X1 <- DATA[[s]]$X
		# X2 <- DATA[[s]]$Y

		# X1.bar <- ( 1 / n1 ) * t(X1) %*% matrix(1,n1,1)

		# H1 <- diag(rep(1,n1)) - (1 / n1) * matrix(1,n1,n1)

		# S1 <- (1/n1) * t(X1) %*% H1 %*% X1

		# X2.bar <- ( 1 / n2 ) * t(X2) %*% matrix(1,n2,1)

		# H2 <- diag(rep(1,n2)) - (1 / n2) * matrix(1,n2,n2)

		# S2 <- (1/n2) * t(X2) %*% H2 %*% X2

		# S <- ( n1*S1 + n2*S2 ) / (n1 + n2 - 2)

		# n <- n1 + n2 - 2

		# S.plus <- ginv(S)

		# D.plus <- t(X1.bar - X2.bar) %*% S.plus %*% (X1.bar - X2.bar)

		# F.plus <- (max(n,p) - min(n,p) + 1) / (n * min(p,n) ) * 1/( 1/n1 + 1/n2 ) * D.plus

		# b.hat <- ( n - 1 ) * ( n + 2) / n^2 * (sum(diag(S))/p)^2 / ( (1/p) * (sum(diag(S %*% S)) - (1/n)*sum(diag(S))^2 ))

		# c.p.n <- sqrt(( (p - n + 1) / ( p + 1) ))

		# Z.plus <- c.p.n * sqrt(n/2) * ( b.hat * F.plus - 1)

		# p.value_z <- 2*( 1 - pnorm(abs(Z.plus)) )

		# p.value_f <- 1 - pf(F.plus,min(n,p),max(n,p) - min(n,p) + 1)
		
		# F.plus.values[s] <- F.plus
		# p.values_f[s] <- p.value_f
		# Z.plus.values[s] <- Z.plus
		# p.values_z[s] <- p.value_z

		# }

# results <- cbind(F.plus.values,p.values_f,Z.plus.values,p.values_z)
# return(results)

# }

# Sriv.onesample.test<-function(X)
# {
	
	# n <- length(X[,1])
	# p <- length(X[1,])

	# X.bar <- ( 1 / n ) * t(X) %*% matrix(1,n,1)

	# H <- diag(rep(1,n)) - (1 / n) * matrix(1,n,n)

	# S <- (1/(n-1)) * t(X) %*% H %*% X

	# S.plus <- ginv(S)

	# D.plus  <- (n-1)*t(X.bar) %*% S.plus %*% X.bar

	# F.plus <- (max(n,p) - min(n,p) + 1) / (n *min(n,p) ) * D.plus

	# b.hat <- ( n - 1 ) * ( n + 2) / n^2 * (sum(diag(S))/p)^2 / ( (1/p) * (sum(diag(S %*% S)) - (1/n)*sum(diag(S))^2 ))

	# #Z.plus <- sqrt(n/2) * ( p * 1/( p - n + 1 ) * b.hat * F.plus - 1)

	# c.p.n <- sqrt(( (p - n + 1) / ( p + 1) ))

	# Z.plus <- c.p.n * sqrt(n/2) * ( b.hat * F.plus - 1)

	# p.value_z <- 2*( 1 - pnorm(abs(Z.plus)))

	# p.value_f <- 1 - pf(F.plus,min(p,n), max(p,n) - min(n,p) + 1)

	# return(c(Z.plus,p.value_z,F.plus,p.value_f))

# }


###############################################################################################
###############################################################################################
################### CHEN AND QIN TEST #########################################################
###############################################################################################
###############################################################################################


Tstat<-function(X1,X2)
{
	
	n1 <- length(X1[,1])
	n2 <- length(X2[,1])
	p <- length(X1[1,])
	
	P1 <- ( sum(X1 %*% t(X1)) - sum(diag( X1 %*% t(X1))) ) / (n1*(n1-1))
	P2 <- ( sum(X2 %*% t(X2)) - sum(diag( X2 %*% t(X2))) ) / (n2*(n2-1))
	P3 <- -2 * sum(X1 %*% t(X2)) / (n1*n2)
	
	T <- P1 + P2 + P3 
	return(T)
}


CVind<-function(n){

	l<-1
	CVind<-matrix(1,n*(n-1)/2,n)
	J<-K<-numeric(n*(n-1)/2)
    for(j in 1:(n-1))
    {
        for(k in (j+1):n)
        {
            CVind[l,c(j,k)]<-0
            J[l]<-j
            K[l]<-k
            l <- l + 1
        }	
    }
	return(list(CVind = CVind , J = J , K = K))
}	


tr.hat<-function(X,CV)
{
	
	accum <- 0
	n <- length(X[,1])
    
	for(l in 1:(n*(n-1)/2))
	{
		xbar_jk <- t( crossprod(CV$CVind[l,], X) / (n-2) )
		a <- crossprod(X[CV$J[l],],X[CV$K[l],])
		b <- crossprod(X[CV$J[l],],xbar_jk)
		c <- crossprod(X[CV$K[l],],xbar_jk)
        
		accum <- accum + (a-b)*(a-c)
	}
	
	return(2 * accum / (n*(n-1)) )
}


CVind12<-function(n1,n2)
{
	l<-1
	J<-numeric(n1)
	K<-numeric(n2)
	CV1<-matrix(1,n1*n2,n1)
	CV2<-matrix(1,n1*n2,n2)
	
	for(j in 1:n1)
	{	
		for(k in 1:n2)
		{
			CV1[l,j]<-0
			CV2[l,k]<-0
			J[l] <- j
			K[l] <- k
			l <- l + 1
		}		
	}
    
	return(list( CV1 = CV1 , CV2 = CV2 , K = K , J = J ))
}


tr12.hat<-function(X1,X2,CV12,n1,n2)
{
	accum <- 0
	
	for(l in 1:(n1*n2))
	{
        
		xbar_j <- t( crossprod(CV12$CV1[l,], X1) / (n1-1) )
		xbar_k <- t( crossprod(CV12$CV2[l,], X2) / (n2-1) )
		a <- crossprod( X1[CV12$J[l],] , X2[CV12$K[l],] )
		b <- crossprod( X1[CV12$J[l],] , xbar_k )
		c <- crossprod( X2[CV12$K[l],] , xbar_j )
        
		accum <- accum + (a - b) * (a - c)
        
	}
	
	return( accum / (n1*n2) )
}

sigma.n <-function(X1,X2,CV1,CV2,CV12,n1,n2,p)
{

	tr12 <- tr12.hat(X1,X2,CV12,n1,n2)
	
	X1bar <- apply(X1,2,mean)
	X2bar <- apply(X2,2,mean)
	
	tr1 <- .C("tracehat",	PACKAGE="highD2pop", as.double(t(X1)),
    as.double(X1bar),
    as.integer(CV1$J-1),
    as.integer(CV1$K-1),
    as.integer(n1),
    as.integer(p),
    result = double(1),
    xbar_jk = double(p))[['result']]
    
	tr2 <- .C("tracehat",	PACKAGE="highD2pop" , as.double(t(X2)),
    as.double(X2bar),
    as.integer(CV2$J-1),
    as.integer(CV2$K-1),
    as.integer(n2),
    as.integer(p),
    result = double(1),
    xbar_jk = double(p))[['result']]
	
	sigma.n <- sqrt( 2 * tr1 / (n1*(n1 - 1)) + 2 * tr2 / (n2*(n2-1)) + 4 * tr12 / (n1*n2) )
	
	return(sigma.n)
}

preChenQin <- function(n1,n2)
{
	CV1<-CVind(n1)
	CV2<-CVind(n2)
	CV12<-CVind12(n1,n2)
	
	return(list(CV1 = CV1, CV2 = CV2, CV12 = CV12))
}


ChenQin.test <- function(X,Y)
{
    n1 <- length(X[,1])
    n2 <- length(Y[,2])
    p <- length(X[1,])
    
    preCQ <- preChenQin(n1,n2)
    
    ChQ <- Tstat(X,Y)/sigma.n(X,Y,preCQ$CV1,preCQ$CV2,preCQ$CV12,n1,n2,p)
    pvalue <- 2*(1 - pnorm(abs(ChQ)))
    return(list(ChQ = ChQ , pvalue = pvalue))
}

ChenQin.sim <- function(DATA)
{
	S <- length(DATA)
	ChQvalues <- pvalues <- numeric(S)
	
	n <- DATA[[1]]$n
	m <- DATA[[1]]$m
	p <- DATA[[1]]$p
	
	preChQ <- preChenQin(n,m)
	
	for(s in 1:S)
    {
        X <- DATA[[s]]$X
        Y <- DATA[[s]]$Y
        
        ChQ <- Tstat(X,Y)/sigma.n(X,Y,preChQ$CV1,preChQ$CV2,preChQ$CV12,n,m,p)
        pvalues[s] <- 2*(1 - pnorm(abs(ChQ)))
        ChQvalues[s] <- ChQ
 
    }
    
	results <- cbind(ChQvalues,pvalues)
	return(results)
}


###############################################################################################
###############################################################################################
################### PCT TEST  #################################################################
###############################################################################################
###############################################################################################

# # PCT.test<- function(X,Y)
# {
	# n <- length(X[,1])
	# m <- length(Y[,1])
	# p <- length(X[1,])	

	
	# X.bar <- ( 1 / n ) * t(X) %*% matrix(1,n,1)
	# Y.bar <- ( 1 / m ) * t(Y) %*% matrix(1,m,1)

	# Sx <- 1 / ( n - 1) * ( t(X)%*%X  - n * X.bar%*%t(X.bar))
	# Sy <- 1 / ( m - 1) * ( t(Y)%*%Y  - m * Y.bar%*%t(Y.bar))
	
	# Sp <- ( ( n - 1 ) * Sx + ( m - 1 ) * Sy) / (n + m - 2)
	
	# Sp.c <- diag(Sp)
	
	# Q <- ( 1 / p ) * sum( (X.bar - Y.bar)^2 / Sp.c ) * ( n * m )/( n + m )
		
		
	# Rp <- diag( 1/sqrt( Sp.c ) ) %*% Sp %*% diag( 1/sqrt( Sp.c ) ) 
	
	# above <- row(Rp) < col(Rp)
	
	# r2 <- Rp[above]^2
	
	# Fstat <- r2 * ( n + m - 2 ) / (1 - r2)
	
	# Sp.shrunk <- Rp.prob <- matrix(0,p,p)
	
	# Rp.prob[above] <- 1 - pf(Fstat, 1, n + m - 4)
	
	# Sp.shrunk[above] <- Sp[above]
	# Sp.shrunk[which(Rp.prob > .05)] <- 0
	
	# Rp.shrunk <- diag( 1/sqrt( Sp.c ) ) %*% Sp.shrunk %*% diag( 1/sqrt( Sp.c ) ) 
	 
		# thing1 <- (2 / p) * ((n + m - 2) / (n + m - 4))^2 * (n + m - 3)/(n + m - 6)
		
		# thing2 <- (2 / p^2) * sum( 1  +  2 * Rp.shrunk[above]^2 ) * ((n + m)/(n + m - 2))^2 
		
		# thing3 <- - (p - 1) / p * ((n + m - 2)/(n + m - 4))^2
	
	# Var.Q <- thing1 + thing2 + thing3
	
	# E.Q <- ( n + m - 2) / (n + m - 4)
	
	# c <- Var.Q / ( 2 * E.Q)
	
	# d <- 2 * E.Q^2 / Var.Q
	
    # pvalue <- 1 - pgamma(Q / c, shape=d/2 ,scale=2)

	# return(list(Q=Q, pvalue=pvalue))
		
# }

# PCT.sim <- function(DATA)
# {
	
	# S <- length(DATA)
	# pvalues <- Qvalues <- numeric(S)
	
	# for(s in 1:S)
    # {
        # X <- DATA[[s]]$X
        # Y <- DATA[[s]]$Y
        # n <- DATA[[s]]$n
        # m <- DATA[[s]]$m
        # p <- DATA[[s]]$p
        
        # X.bar <- ( 1 / n ) * t(X) %*% matrix(1,n,1)
        # Y.bar <- ( 1 / m ) * t(Y) %*% matrix(1,m,1)
        
        # Sx <- 1 / ( n - 1) * ( t(X)%*%X  - n * X.bar%*%t(X.bar))
        # Sy <- 1 / ( m - 1) * ( t(Y)%*%Y  - m * Y.bar%*%t(Y.bar))
		
        # Sp <- ( ( n - 1 ) * Sx + ( m - 1 ) * Sy) / (n + m - 2)
		
        # Sp.c <- diag(Sp)
        
        # Q <- ( 1 / p ) * sum( (X.bar - Y.bar)^2 / Sp.c ) * ( n * m )/( n + m )
		
		
        # Rp <- diag( 1/sqrt( Sp.c ) ) %*% Sp %*% diag( 1/sqrt( Sp.c ) ) 
        
        # above <- row(Rp) < col(Rp)
        
        # r2 <- Rp[above]^2
        
        # Fstat <- r2 * ( n + m - 2 ) / (1 - r2)
        
        # Sp.shrunk <- Rp.prob <- matrix(0,p,p)
        
        # Rp.prob[above] <- 1 - pf(Fstat, 1, n + m - 4)
        
        # Sp.shrunk[above] <- Sp[above]
        # Sp.shrunk[which(Rp.prob > .05)] <- 0
        
        # Rp.shrunk <- diag( 1/sqrt( Sp.c ) ) %*% Sp.shrunk %*% diag( 1/sqrt( Sp.c ) ) 
        
        # thing1 <- (2 / p) * ((n + m - 2) / (n + m - 4))^2 * (n + m - 3)/(n + m - 6)
		
        # thing2 <- (2 / p^2) * sum( 1  +  2 * Rp.shrunk[above]^2 ) * ((n + m)/(n + m - 2))^2 
		
        # thing3 <- - (p - 1) / p * ((n + m - 2)/(n + m - 4))^2
        
        # Var.Q <- thing1 + thing2 + thing3
        
        # E.Q <- ( n + m - 2) / (n + m - 4)
        
        # c <- Var.Q / ( 2 * E.Q)
        
        # d <- 2 * E.Q^2 / Var.Q
        
        # pvalues[s] <- pgamma(Q / c, shape=d/2 ,scale=2)
        
        # Qvalues[s] <- Q
        
    # }
    
	# results <- cbind(Qvalues,pvalues)
	# return(results)
# }

