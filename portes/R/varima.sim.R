"varima.sim" <-
function (phi = NULL, theta = NULL, d = NA, sigma, n, constant = NA, 
         trend = NA, demean = NA, innov = NULL, innov.dist=c("normal","t","stable"),  
         df=1, StableParameters=NA, Trunc.Series = NA) 
{
	   innov.dist <- match.arg(innov.dist)
	    if (!is.null(phi) && class(phi)!="array" && class(phi)!="numeric")
		stop("Phi must be enterd as NULL or array with dimension (k*k*p) or numeric")
	    if (!is.null(theta) && class(theta)!="array" && class(theta)!="numeric")
		stop("Theta must be enterd as NULL or array with dimension (k*k*q) or numeric")
	    sigma <- as.matrix(sigma)
	    k <- NCOL(sigma)
	    if (all(is.na(d))||all(d==0)) 
	        d <- rep(0, k)
	    if (length(d) != k) 
	      stop("d must be entered as a vector with length equal to number of sigma rows")
 	    if (any(d < 0)) 
            stop("number of differences must be a nonnegative integer/integers")
	    if (all(is.na(constant))) 
 	      constant <- rep(0, k)
	    if (length(constant) != k) 
	      stop("constant must be entered as a vector with length equal to number of sigma rows")
	    if (all(is.na(trend))) 
      	  trend <- rep(0, k)
	    if (length(trend) != k) 
	        stop("trend must be entered as a vector with length equal to number of sigma rows")
	    if (all(is.na(demean))) 
	        demean <- rep(0, k)
	    if (length(demean) != k) 
	        stop("demean must be entered as a vector with length equal to number of sigma rows")
	    if (class(phi) == "numeric")
	        phi <- array(phi,dim=c(1,1,length(phi)))
	    if (class(theta)=="numeric")
	        theta <- array(theta,dim=c(1,1,length(theta)))
 	    if (all(phi == 0))
 	        phi <- NULL
	    if (all(theta == 0))
              theta <- NULL
          p <- ifelse(is.null(phi),0,dim(phi)[3])
          q <- ifelse(is.null(theta),0,dim(theta)[3])
          if (p > 0 && ((dim(phi)[1] != dim(phi)[2])||dim(phi)[2] != k))
             stop("Wrong dimensions of phi or/and sigma")
          if (q > 0 && ((dim(theta)[1] != dim(theta)[2])||dim(theta)[2] != k))
             stop("Wrong dimensions of theta or/and sigma")
          if(is.null(innov)&& innov.dist == "stable"){
              StableQ <- all(!is.na(StableParameters))
              if (StableQ) {
                  StableParameters <- matrix(StableParameters,nrow=k)
               if (NCOL(StableParameters) != 4) 
                  stop("StableParameters must be a numeric vector/matrix with 4 stable parameters") 
               ALPHA<-StableParameters[,1]
               BETA<-StableParameters[,2]
               GAMMA<-StableParameters[,3]
               DELTA<-StableParameters[,4]
              }
             else 
               stop("StableParameters must be a numeric vector/matrix with 4 stable parameters")
           }
          if (p == 0) {
             if(!is.null(innov)){ 
                 innov <- as.matrix(innov)
                 epsilon <- innov
                 stopifnot (NROW(epsilon) == n && NCOL(epsilon) == k)
             }
             else if(is.null(innov)&& innov.dist == "normal")
                  epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))
             else if(is.null(innov)&& innov.dist == "t")
                  epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(n + q)),ncol=n + q)))/sqrt(rchisq(k*(n + q),df)/df)
             else if(is.null(innov)&& innov.dist == "stable")
                  epsilon <- rStable(n + q, ALPHA, BETA, GAMMA, DELTA)
             if (q == 0){  ## Simulate white noise
                    SlopDrift <- t(matrix(trend * rep(1:NROW(epsilon),each=k),nrow=k,ncol=NROW(epsilon)))
                 DriftTerm <- sweep(SlopDrift, 2L, -constant, check.margin = FALSE)
                 CenterData <- scale(epsilon, center = -demean, scale = FALSE)
                 Sim.Series <- DriftTerm + CenterData
                 for (i in 1:length(d)){
                   if(d[i]>0)
                     Sim.Series[,i] <- as.matrix(diffinv(Sim.Series[,i], differences = d[i]))[-(1:d[i]),]
                   else
                     Sim.Series[,i] <- Sim.Series[,i]   
                   }
                 return(ts(Sim.Series))
               }
              else
                 InvertQ(theta) ## Simulate VMA(q)
                 psi <- array(c(diag(k), theta), dim = c(k, k, q + 1))
                 if(!is.null(innov))
                    epsilon <- ts(rbind(innov,as.matrix(innov[sample(x=1:n,size=q,replace=TRUE),])))
                 Sim.VMA <- vma.sim(psi = psi, a = epsilon)
                 SlopDrift <- t(matrix(trend * rep(1:NROW(Sim.VMA),each=k),nrow=k,ncol=NROW(Sim.VMA)))
                 DriftTerm <- sweep(SlopDrift, 2L, -constant, check.margin = FALSE)
                 CenterData <- scale(Sim.VMA, center = -demean, scale = FALSE)
                 Sim.Series <- DriftTerm + CenterData
                 for (i in 1:length(d)){
                    if(d[i]>0)
                      Sim.Series[,i] <- as.matrix(diffinv(Sim.Series[,i], differences = d[i]))[-(1:d[i]),]
                    else
                      Sim.Series[,i] <- Sim.Series[,i]   
                 }
                return(ts(Sim.Series))
              }
          else
            if (is.na(Trunc.Series)) 
               Trunc.Series <- min(100,ceiling(n/3))
            FirstSim.Series <- matrix(numeric(0), nrow = n, ncol = k)
            r <- max(p, q)
            psi <- ImpulseVMA(phi = phi, theta = theta, Trunc.Series = Trunc.Series)
            if(!is.null(innov)) 
                epsilon <- matrix(innov[sample(x=1:n,size=Trunc.Series + r,replace=TRUE),],nrow=Trunc.Series + r,ncol=k)             
            else if(is.null(innov)&& innov.dist == "normal")
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(Trunc.Series + r)),ncol=Trunc.Series + r)))
            else if(is.null(innov)&& innov.dist == "t")
                epsilon <- t(crossprod(chol(sigma),matrix(stats::rnorm(k*(Trunc.Series + r)),ncol=Trunc.Series + r)))/sqrt(rchisq(k*(Trunc.Series + r),df)/df)
            else if(is.null(innov)&& innov.dist == "stable")
                epsilon <- rStable(Trunc.Series + r, ALPHA, BETA, GAMMA, DELTA)
            FirstSim.Series[1:r, ] <- vma.sim(psi = psi, a = epsilon)
            a <- matrix(epsilon[1:r, ], nrow = r, ncol = k)
            if(!is.null(innov)) 
               epsilon <- rbind(a,innov)            
            else if(is.null(innov)&& innov.dist == "normal")
               epsilon <- rbind(a,t(crossprod(chol(sigma),matrix(stats::rnorm(k*n),ncol=n))))
            else if(is.null(innov)&& innov.dist == "t")
               epsilon <- rbind(a,t(crossprod(chol(sigma),matrix(stats::rnorm(k*n),ncol=n)))/sqrt(rchisq(k*n,df)/df))
            else if(is.null(innov)&& innov.dist == "stable")
               epsilon <- rbind(a,rStable(n, ALPHA, BETA, GAMMA, DELTA))
            if (q > 0) { ## Simulate VARMA(p,q)
                extend.psi <- array(c(diag(k), theta, rep(0, k * k *(n - q))), dim = c(k, k, n + 1))
                u <- matrix(numeric(0), nrow = n, ncol = k)
                for (i in (q + 1):(n + q)) {
                   out <- 0
                   for (j in 0:q) {
                      out = out + crossprod(t(extend.psi[, , j + 1]),epsilon[i - j, ])
                   }
                  u[i - q, ] <- out
                }
             }
            else 
               u <- epsilon  ## Simulate VAR(p)
               for (i in (r + 1):n) {
                 temp2 <- 0
                 for (j in 1:p) temp2 <- temp2 + crossprod(t(phi[, , j]),FirstSim.Series[i - j, ])
                   FirstSim.Series[i, ] <- temp2 + u[i, ]
               }
              SlopDrift <- t(matrix(trend * rep(1:NROW(FirstSim.Series),each=k),nrow=k,ncol=NROW(FirstSim.Series)))
              DriftTerm <- sweep(SlopDrift, 2L, -constant, check.margin = FALSE)
              CenterData <- scale(FirstSim.Series, center = -demean, scale = FALSE)
              Sim.Series <- DriftTerm + CenterData
            for (i in 1:length(d)){
              if(d[i]>0)
                 Sim.Series[,i] <- as.matrix(diffinv(Sim.Series[,i], differences = d[i]))[-(1:d[i]),]
              else
                 Sim.Series[,i] <- Sim.Series[,i]   
             }
          return(ts(Sim.Series))
}




