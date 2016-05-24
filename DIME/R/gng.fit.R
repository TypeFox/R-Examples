gng.fit <-
function(data, avg = NULL, K = 2, weights = NULL, weights.cutoff = -1.345,
  pi = NULL, mu = NULL, sigma = NULL, beta = NULL, tol=1e-5, 
  max.iter=2000, th= NULL)
{
  x <- unlist(data);
  n <- length(x);		
	if(is.null(weights)) weights <- rep(1,length(x));
  if(!is.null(weights) & !is.character(weights) & length(weights)!=length(x)){
      return(cat('Error: please input weights need to be the same length
       as data\n'));
  }
  if(!is.null(weights) && is.character(weights) && is.null(avg)){
      return(cat('Error: When using weights, please input the mean 
        (or log intensities) of data\n'));
  }
  if(!is.null(weights) && is.character(weights)){
    weights <- match.arg(tolower(weights),c("lower","upper","full"));
	  weights <- switch(weights,
		lower = huber(avg, weights.cutoff,'lower'),
		upper = huber(avg, weights.cutoff,'upper'),
		full = huber(avg, weights.cutoff,'full'))
  }
        
	if(K < 1) stop("It is expected that there is at least one normal component; therefore, K needs to be greater than 0.");
	# let C number of component to be C = K + 2
	C <- K + 2 ;       
	
	# initialize parameters for the exponential components
	if(is.null(th)) {
		if (length(x[x<0])>0){
		if(abs(max(x[x<0]))< 2){
		   th1<-runif(1,abs(max(x[x<0])),2);
		   }else{
       th1<-runif(1,0,2);
       }
    }else{
    th1 <- runif(1,0,2);
    }
		if (length(x[x>0]) > 0){
    if(abs(min(x[x>0]))< 2){
       th2<-runif(1,abs(min(x[x>0])),2);
    }else{
       th2<-runif(1,0,2);
		}
		}else{
		 th2<-runif(1,0,2);
	 }
 }else { 
		th1<-abs(max(-th,min(x)));
		th2<-min(th,max(x));
 	}	
	if(is.null(beta)){                
		# initialize value for beta estimate for negatives
		beta <-c(0,0);
		if (length(x[x<0])>0){
    meanX_neg <- mean(x[x<0]);
		stdevX_neg<- sd(x[x<0]);
           	beta[1] <- abs(rnorm(1,meanX_neg,stdevX_neg));
           	
    }else{
         beta[1] <- runif(1,20,30);
    }
    if (length(x[x>0])>0){
		# initialize value for beta estimate for positives
		meanX_pos <- mean(x[x>0]);
		stdevX_pos <- sd(x[x>0]);
    beta[2] <- abs(rnorm(1,meanX_pos,stdevX_pos));
		}else{ 
    beta[2] <- abs(rnorm(1,mean(x),sd(x)));
    }
	}

	# initialize the proportion of mixture components
	if(is.null(pi)){
      	pi <- runif(C,0,1);
            pi <- pi/sum(pi);        
	}        

	# initialize parameters for normal components
	if(is.null(mu)){
		mu <- sample(seq(-1,1,length = 100), K, replace=TRUE);	
        }        
	if(is.null(sigma)){
	      sigma <-rep(1/(K),K);
	} 
 	# get indicator matrix for negative and positive x's
      I1 <- (x < (-th1))+0;        
	I2 <- (x > th2)+0;
ENK <- .C("ENKE", input = as.double(x), length = as.integer(n), components = as.integer(C), weights = as.double(weights),
	 max.iterations = as.integer(max.iter), zero = as.double(tol), pi = as.double(pi), mu = as.double(mu), 
	 sigma = as.double(sigma), beta = as.double(beta), negatives = as.integer(I1), positives = as.integer(I2), 
	 thresholds = as.double(c(th1, th2)), iterations = integer(1),PACKAGE="DIME");
if (!is.nan(ENK$beta[1]) && !is.nan(ENK$beta[2])){
phi <- matrix(0, n, C);
for (k in 2:(C-1)){
	phi[,k]<-ENK$pi[k]*dnorm(x, ENK$mu[k-1], ENK$sigma[k-1]);
	}
	
phi[,1] <- ENK$pi[1]*(I1)*dexp((-x-th1),rate = 1/ENK$beta[1]);
phi[,C] <- ENK$pi[C]*(I2)*dexp((x-th2),rate = 1/ENK$beta[2]);

# summing all columns of phi
sum1 <- c(phi %*% matrix(rep(1,ncol(phi))));

loglike <- sum(weights * log(sum1));
AIC <- loglike - (3*K+3);
BIC <- 2*loglike - (3*K+3)*log(n);

# calculating the estimated model for undifferentiated probes
if (C > 3) { # no rowsums is needed for C<=3
	f0_psi <- rowSums(phi[,(2):(C-1)]);
}
else {
	f0_psi <- phi[,2];
}

# calculating the estimated model for all probes
f <- rowSums(phi[,1:C]);
# calculating the local fdr
if (C > 3) { # no rowsums is needed for C<=3
	fdr <- f0_psi/(f*sum((ENK$pi[2:(C-1)])));
}
else {
	fdr <-f0_psi/(f*ENK$pi[2]);
}

param_estim <- list(fdr=fdr, beta=ENK$beta,mu=ENK$mu,sigma=ENK$sigma,
    loglike=loglike,AIC=AIC, BIC=BIC, iter=ENK$iterations,
    range=c(min(x),max(x)),th1=th1,th2=th2, K = K,pi=ENK$pi,phi=phi,name="GNG");
}else{ # if is.nan(ENK$beta)
AIC = -Inf;
BIC = -Inf;
param_estim <- list(AIC = AIC,BIC = BIC);
}
return(param_estim);
}

