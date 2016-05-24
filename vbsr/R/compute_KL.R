compute_KL = function(Zmat,alpha,pval){
  #define function to sample from truncated normal distribution.
  rtnorm <- function(n,mu,sigma,a,b){
    accept <- rep(0,n);
    x <- rep(0,n);
    u <- rep(0,n);
    if(a>mu){
      l<-a;
    } else if(b<mu){
      l<-b;
    } else{
      l<-mu;
    }
    wh <- 1:n;
    while(sum(accept)<n){
      x[wh] <- runif(length(wh),a,b);
      u[wh] <- runif(length(wh),0,1);
      ntrue <- which(u[wh]<dnorm(x[wh],mu,sigma)/dnorm(l,mu,sigma));
      accept[wh][ntrue]<-TRUE;
      wh<-which(accept==FALSE);
    }
    return(x);
  }
  
  
  
	n <- nrow(Zmat);
	m <- ncol(Zmat);
	kl_vec <- rep(0,m);
	kl_vec2 <- rep(0,1e3);

	if(pval <1){
		m2 <- n*(1/pval);
		cc <- -qnorm(1/(2*m2));
		a <- -qnorm(pval/2);
		b <- -qnorm(alpha/2);
		sig_ref <- 1 + (a*dnorm(a)-b*dnorm(b))/(pnorm(b)-pnorm(a));
		#print(c(a,b));
		for (i in 1:m){

			wvec <- intersect(which(abs(Zmat[,i])<b),which(abs(Zmat[,i])>a));
			#print(length(wvec));
			#wvec <- 1:m;
			mu <- mean(Zmat[wvec,i]);
			sig <- var(Zmat[wvec,i]);
			kl_vec[i] <- mu^2/2 + 0.5*(sig/sig_ref-1-log(sig/sig_ref));
		}

		for (i in 1:1e3){
			vec3 <- rep(0,n);
			vec4 <- rep(0,n);
			#wvec2 <- intersect(which(abs(vec3)<b),which(abs(vec3)>a));
			#wvec2 <- 1:m;
			z <- rbinom(n,1,.5);
			vec3 <- rtnorm(n,0,1,a,cc);
			vec4 <- rtnorm(n,0,1,-cc,-a);
			vec3[which(z==1)]<-vec4[which(z==1)];
			wvec2 <- intersect(which(abs(vec3)<b),which(abs(vec3)>a));
			mu2 <- mean(vec3[wvec2]);
			sig2 <- var(vec3[wvec2]);
			kl_vec2[i] <- mu2^2/2 + 0.5*(sig2/sig_ref-1-log(sig2/sig_ref));
		}

	}else {
		a <- -qnorm(alpha/2);
		sig_ref <- (1-2*a*dnorm(a))/(pnorm(a)-pnorm(-a));
		for (i in 1:m){
			wvec <- which(abs(Zmat[,i])<a);
			mu <- mean(Zmat[wvec,i]);
			sig <- var(Zmat[wvec,i]);
			kl_vec[i] <- mu^2/2 + 0.5*(sig/sig_ref-1-log(sig/sig_ref));
		}
		for (i in 1:1e3){
			vec3 <- rnorm(n);
			wvec2 <- which(abs(vec3)<a);
			#wvec2 <- 1:m;
			mu2 <- mean(vec3[wvec2]);
			sig2 <- var(vec3[wvec2]);
			kl_vec2[i] <- mu2^2/2 + 0.5*(sig2/sig_ref-1-log(sig2/sig_ref));
		}

	}

	res_list <- vector("list",4);
	names(res_list) <- c("kl_vec","min_kl","mean_kl","se_kl");

	res_list$kl_vec <- log(kl_vec);
	res_list$min_kl <- min(log(kl_vec));
	res_list$mean_kl <- mean(log(kl_vec2));
	res_list$se_kl <- sd(log(kl_vec2));

	return(res_list);

}
