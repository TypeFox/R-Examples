vbsr = function(y,
		X,
		ordering_mat=NULL,
		eps=1e-6,
		exclude=NULL,
		add.intercept=TRUE,
		maxit = 1e4,
		n_orderings = 10,
    family = "normal",
		scaling = TRUE,
		return_kl = TRUE,
		estimation_type = "BMA",
		bma_approximation = TRUE,
		screen = 1.0,
		post=0.95,
		already_screened = 1.0,
		kl = 0.99,
		l0_path=NULL,
    cleanSolution=FALSE){

	n <- nrow(X);
	m <- ncol(X);
	if(add.intercept==TRUE){
		X <- cbind(rep(1,n),X);
		m <- m+1;
	}
	if(!is.null(post)){
		path_length=1;
		l0_path=-(qchisq(0.05/m,1,lower.tail=FALSE)-log(n)+2*log((1-post)/(post)));
	}else{
    path_length=length(l0_path)
    if(path_length==0){
      stop("invalid penalty parameter path specification")
    }
	}
	#n <- nrow(X);

	#m <- ncol(X);
	#print(c("m: ",m));
	#eps = 1e-6;
	#path_length = 50;
	if(is.null(exclude)){
		if(add.intercept==TRUE){
			exclude = rep(0,m);
			exclude[1] = 1;
		}else{
			exclude = rep(0,m);
		}
	}else {
		if(add.intercept==TRUE){
			exclude <- c(1,exclude);
			if(length(exclude)!=m){
				stop("Non-penalization indicator vector of wrong length: ",length(exclude)-1,"!=",m-1);
			}
			if(length(which(exclude==0))+length(which(exclude==1))<m){
				stop("Not all elements of non-penalization indicator vector are 1 or 0.");
			}

		} else {
			if(length(exclude)!=m){
				stop("Non-penalization indicator vector of wrong length: ",length(exclude),"!=",m);
			}
			if(length(which(exclude==0))+length(which(exclude==1))<m){
				stop("Not all elements of non-penalization indicator vector are 1 or 0.");
			}

		}
	}



	#maxit = 1e4;
	#n_orderings = 10;
	if(family == "normal"){
		regress <- 1;
	} else if (family == "binomial"){
		regress <- 0;
	} else {
		stop("Improper type of regression provided. Must be either 'normal' or 'binomial'.");
	}
	#regress = 0;

	if(scaling==TRUE){
		scale <- 1;
	} else if (scaling==FALSE){
		scale <- 0;
	} else {
		stop("Improper design matrix scaling parameter provided.");
	}

	#scale = 1;

	if(return_kl==TRUE){
		error <- 1;
	} else if (return_kl==FALSE){
		error <- 0;
	} else {
		stop("Improper KL switch for returning KL results.");
	}
	#error = 1; 


	if(estimation_type=="BMA"){
		est <- 1;
	} else if (estimation_type =="MAXIMAL"){
		est <- 0;
	} else {
		stop("Improper global estimation type.  Must be either 'BMA' or 'MAXIMAL'.");
	}

	#est = 1;

	if(bma_approximation==TRUE){
		approx <- 1;
	} else if (bma_approximation==FALSE){
		approx <- 0;
	} else {
		stop("Improper Bayesian model averaging z-score approximate estimation indicator.");
	}
	#approx = 1;

	
	if(screen > 1 || screen <= 0){
		stop("Improper marginal screening parameter.");
	}
	
	if(already_screened > 1 || already_screened <= 0){
		stop("Improper previously pre-screened parameter.");
	}
	#screen = .1;
	#already_screened = 1;


	#kl = 0.99;

	if(kl >= 1 || kl <=0){
		stop("Improper KL percentile.");
	}

	total_replicates = path_length*n_orderings;
	var_y = var(y);

	beta_chi_mat = double(m*path_length);
	beta_mu_mat = double(m*path_length);
	beta_sigma_mat = double(m*path_length);
	e_beta_mat = double(m*path_length);
	beta_p_mat = double(m*path_length);
	lb_mat = double(n_orderings*path_length);
	kl_mat = double(n_orderings*path_length);
	#print(c(n,m));
	
	beta_chi = double(m);
	beta_mu = double(m);
	beta_sigma = double(m);
	beta_p = double(m);
	lb1 = double(1);
	
	sma<-.C("run_marg_analysis",
			as.double(eps),
			as.integer(exclude),
			as.integer(maxit),
			as.integer(regress),
			as.integer(scale),
			as.double(X),
			as.double(y),
			as.double(var_y),
			as.integer(n),
			as.integer(m),
			beta_chi,
			beta_mu,
			beta_sigma,
			beta_p,
			lb1,PACKAGE="vbsr");
	
	wexc <- which(exclude==1);

	#compute score vector for adaptively reweighted logistic function
	score_vec <- sort(sma[[11]][-wexc]^2+log(sma[[13]][-wexc]),decreasing=T);
	#print(score_vec);
	#print(c("MAX:",max(score_vec,33)));
	#determine automatic path based on data
	if(is.null(l0_path)){
		if(sqrt(n)>m){
			#print('l0 larger');
			l0_path = seq(-score_vec[1],-min(score_vec),length.out=path_length);
		}else{
			#print('l0 smaller');
			l0_path = seq(-score_vec[1],-score_vec[round(sqrt(n))],length.out=path_length);
		}
		#print(l0_path);
	}
	l0_path <- rev(l0_path);
	pb_path = 1/(1+exp(-0.5*l0_path));
	#p_est <- rep(0,n_orderings);

	#if(max_l0==TRUE){
	#	l0_path <- 0;
	#	pb_path <- 1/(1+exp(-0.5*l0_path));
	#	path_length <- 1;
	#	max_pb <- 1;
	#	#if(path_length==1){error <- 0;}
	#}else{
	#	max_pb <- 0;
	#}
	#which always keep
	#need to fix kl statistic function to accept path lengths of 1
	if(path_length==1){error <- 0;}	

	#compute sma p-values if pre-screening:
	if(screen <1){
		sma_p <- exp(pchisq(sma[[11]]^2,1,log.p=T,lower.tail=F));
		#print(sma_p);
		wkeep <- which(sma_p<screen);
		wkeep <- sort(unique(c(wexc,wkeep)));
		#print(wkeep);
		if(length(wkeep)==length(wexc)){
			stop("All features removed during marginal screening, raise marginal screening threshold.\n");
		}
		m_s <- length(wkeep);
		dv <- 1:m_s;
		names(dv) <- as.character(wkeep);
		wexc_s <- as.numeric(as.character(dv[wexc]));
		exclude_s <- rep(0,m_s);
		exclude_s[wexc_s] <- 1;
		#print(wexc_s);
		X <- X[,wkeep];
		m <- m_s;
		beta_chi_mat = double(m*path_length);
		beta_mu_mat = double(m*path_length);
		beta_sigma_mat = double(m*path_length);
		e_beta_mat = double(m*path_length);
		beta_p_mat = double(m*path_length);
		exclude <- exclude_s;
		#already_screened<-screen;
		#print(dim(X));
		#print(m);

		penalty_factor = rep(1,m);
		if(is.null(ordering_mat)){
			ordering_mat = matrix(0,m,n_orderings);
			for(i in 1:n_orderings){
				ordering_mat[,i]<- sample(1:m,m);
			}
			ordering_mat <- ordering_mat -1;
		}

		if(!is.null(ordering_mat)){
			if(nrow(ordering_mat)!=m || ncol(ordering_mat)!= n_orderings){
				stop("Provided orderings are of wrong dimension.\n");
			}
			if(min(ordering_mat)>0){
				stop("Index must start from 0.\n");
			}
			for(i in 1:n_orderings){
				if(sd(sort(ordering_mat[,i])-(0:(m-1)))>0){
					stop("Incorrect ordering:",i);
				}
			}
		}
		#gc();
	} else {
		penalty_factor = rep(1,m);
		if(is.null(ordering_mat)){
			ordering_mat = matrix(0,m,n_orderings);
			for(i in 1:n_orderings){
				ordering_mat[,i]<- sample(1:m,m);
			}
			ordering_mat <- ordering_mat -1;
		}
		if(!is.null(ordering_mat)){
			if(nrow(ordering_mat)!=m || ncol(ordering_mat)!= n_orderings){
				stop("Provided orderings are of wrong dimension.\n");
			}
			if(min(ordering_mat)>0){
				stop("Index must start from 0.\n");
			}
			for(i in 1:n_orderings){
				if(sd(sort(ordering_mat[,i])-(0:(m-1)))>0){
					stop("Incorrect ordering:",i);
				}
			}
		}
	}
	#if(!is.null(n_threads)){
	#  nthreads=n_threads;
	#}else{
	  nthreads=1;
	#}
	result <- c();
	while(length(result)==0){
		try(result<-.C("run_vbsr_wrapper",
			as.double(eps),
			as.double(l0_path),
			as.double(pb_path),
			as.integer(exclude),
			as.double(penalty_factor),
			as.integer(maxit),
			as.integer(path_length),
			as.integer(n_orderings),
			as.integer(regress),
			as.integer(scale),
			as.integer(est),
			as.integer(error),
			as.double(kl),
			as.integer(approx),
			as.integer(total_replicates),
			as.double(X),
			as.double(y),
			as.double(var_y),
			as.integer(n),
			as.integer(m),
			as.integer(ordering_mat),
			as.double(beta_chi_mat),
			as.double(beta_mu_mat),
			as.double(beta_sigma_mat),
			as.double(e_beta_mat),
			as.double(beta_p_mat),
			as.double(lb_mat),
			as.double(kl_mat),
			as.integer(nthreads),
			PACKAGE="vbsr"),silent=TRUE);
		if(length(result)==0&&path_length>1){
			#rm(result);
			#gc();
			#result <- c();
			path_length <- path_length-1;
			l0_path <- l0_path[-1];
			pb_path <- pb_path[-1];
			beta_chi_mat = double(m*path_length);
			beta_mu_mat = double(m*path_length);
			beta_sigma_mat = double(m*path_length);
			e_beta_mat = double(m*path_length);
			beta_p_mat = double(m*path_length);
			lb_mat = double(n_orderings*path_length);
			kl_mat = double(n_orderings*path_length);
		} else if (length(result)==0&&path_length<=1){
			stop("solution does not exist for any of path specified");
		}
	}

	if(scale==1&&add.intercept==TRUE){
		mult <- c(1,apply(X[,-1],2,sd)*sqrt((n-1)/n));
	} else if (scale==1){
		mult <- apply(X[,-1],2,sd)*sqrt((n-1)/n);
	}else{
		mult <- rep(1,m);
	}

	if (screen==1){


						

		result_list = vector("list",17);
		names(result_list)=c("beta_chi",
					"beta_mu",
					"beta_sigma",
					"e_beta",
					"beta_p",
					"lb",
					"kl",
					"l0_path",
					"sma_beta",
					"sma_chi",
					"kl_min",
					"kl_mean",
					"kl_se",
					"ordering_mat",
					"which_excluded",
					"kl_index");
		result_list$kl = matrix(result[[28]],n_orderings,path_length);
		result_list$lb = matrix(result[[27]],n_orderings,path_length);
		result_list$beta_p = matrix(result[[26]],m,path_length);
		result_list$e_beta = matrix(result[[25]],m,path_length)/(mult);
		#result_list$e_beta = matrix(result[[25]],m,path_length);
		result_list$beta_sigma = matrix(result[[24]],m,path_length)/(mult^2);
		#result_list$beta_sigma = matrix(result[[24]],m,path_length);
		result_list$beta_mu = matrix(result[[23]],m,path_length)/(mult);
		#result_list$beta_mu = matrix(result[[23]],m,path_length);
		result_list$beta_chi = matrix(result[[22]],m,path_length);
		result_list$l0_path = l0_path;
		result_list$sma_beta = sma[[12]];
		result_list$sma_chi = sma[[11]];
		result_list$ordering_mat = ordering_mat;
		result_list$which_excluded = wexc;

		#compute KL statistics
		if(error == 1){
			kl_res <- compute_KL(result_list$beta_chi[-wexc,],1-kl,already_screened);
			#print(kl_res);
			result_list$kl <- rev(kl_res$kl_vec);
			result_list$kl_min <- kl_res$min_kl;
			result_list$kl_mean <- kl_res$mean_kl;
			result_list$kl_se <- kl_res$se_kl;

			nl <- length(result_list$kl);
			a_list <- vector("list",6);
			c_list <- rep(0,6);
			A <- cbind(result_list$kl[1:(nl-1)],result_list$kl[2:(nl)])
			c_list[1] <- result_list$kl_min;
			c_list[2] <- result_list$kl_mean;
			c_list[3] <- result_list$kl_min+result_list$kl_se;
			c_list[4] <- result_list$kl_mean+result_list$kl_se;
			c_list[5] <- result_list$kl_min+2*result_list$kl_se;
			c_list[6] <- result_list$kl_mean+2*result_list$kl_se;
			b_list <- rep(0,6);
			for(i in 1:6){
				a_list[[i]]<-intersect(which(A[,1]<=c_list[i]),which(c_list[i]<=A[,2]));
				if(length(a_list[[i]])>0){
					mA <- max(a_list[[i]]);
					wm <- which.min(abs(A[mA,]-c_list[[i]]));
					if(wm==1){
						b_list[i] <- mA;
					} else{
						b_list[i] <- mA+1;
					}
				}else{
					if(c_list[i]<=result_list$kl[1]){
						b_list[i] <- 1;
					}else{
						b_list[i] <- length(result_list$kl);
					}
				}
			}
			b_list <- path_length+1-b_list;
			names(b_list)<-c("kl_min",
					"kl_mean",
					"kl_min_1se",
					"kl_mean_1se",
					"kl_min_2se",
					"kl_mean_2se");
			result_list$kl_index <- b_list;

		}

	} else {


		#compute KL statistics


		result_list = vector("list",18);
		names(result_list)=c("beta_chi",
					"beta_mu",
					"beta_sigma",
					"e_beta",
					"beta_p",
					"lb",
					"kl",
					"l0_path",
					"sma_beta",
					"sma_chi",
					"screened_ind",
					"kl_min",
					"kl_mean",
					"kl_se",
					"ordering_mat",
					"which_excluded",
					"kl_index");
		result_list$kl = matrix(result[[28]],n_orderings,path_length);
		result_list$lb = matrix(result[[27]],n_orderings,path_length);
		result_list$beta_p = matrix(result[[26]],m,path_length);
		result_list$e_beta = matrix(result[[25]],m,path_length)/(mult);
		result_list$beta_sigma = matrix(result[[24]],m,path_length)/(mult^2);
		result_list$beta_mu = matrix(result[[23]],m,path_length)/(mult);
		result_list$beta_chi = matrix(result[[22]],m,path_length);
		result_list$l0_path = l0_path;
		result_list$sma_beta = sma[[12]];
		result_list$sma_chi = sma[[11]];
		result_list$screened_ind = wkeep;
		result_list$ordering_mat = ordering_mat;
		result_list$which_excluded = wexc;
		if(error == 1){
			kl_res <- compute_KL(result_list$beta_chi[-wexc,],1-kl,screen);
			
			#print(kl_res);
			result_list$kl <- rev(kl_res$kl_vec);
			result_list$kl_min <- kl_res$min_kl;
			result_list$kl_mean <- kl_res$mean_kl;
			result_list$kl_se <- kl_res$se_kl;
			nl <- length(result_list$kl);
			a_list <- vector("list",6);
			c_list <- rep(0,6);
			A <- cbind(result_list$kl[1:(nl-1)],result_list$kl[2:(nl)])
			c_list[1] <- result_list$kl_min;
			c_list[2] <- result_list$kl_mean;
			c_list[3] <- result_list$kl_min+result_list$kl_se;
			c_list[4] <- result_list$kl_mean+result_list$kl_se;
			c_list[5] <- result_list$kl_min+2*result_list$kl_se;
			c_list[6] <- result_list$kl_mean+2*result_list$kl_se;
			b_list <- rep(0,6);
			for(i in 1:6){
				a_list[[i]]<-intersect(which(A[,1]<=c_list[i]),which(c_list[i]<=A[,2]));
				if(length(a_list[[i]])>0){
					mA <- max(a_list[[i]]);
					wm <- which.min(abs(A[mA,]-c_list[[i]]));
					if(wm==1){
						b_list[i] <- mA;
					} else{
						b_list[i] <- mA+1;
					}
				}else{
					if(c_list[i]<=result_list$kl[1]){
						b_list[i] <- 1;
					}else{
						b_list[i] <- length(result_list$kl);
					}
				}
			}
			b_list <- path_length+1-b_list;
			names(b_list)<-c("kl_min",
					"kl_mean",
					"kl_min_1se",
					"kl_mean_1se",
					"kl_min_2se",
					"kl_mean_2se");
			result_list$kl_index <- b_list;

		}
	}	
  modUnique <- which(!duplicated(round(result_list$lb,-log10(eps)-1)));
  lb1 <- result_list$lb[modUnique];
  #print(result_list$lb);
  modProb <- exp(lb1-max(lb1))/(sum(exp(lb1-max(lb1))));
  
  
	if(!is.null(post)){
    result_list2 <- list();
		result_list2$beta <- result_list$e_beta[-wexc];
    result_list2$alpha <- result_list$e_beta[wexc];
    #result_list2$betaSE <- sqrt(result_list$e_beta[-wexc]^2+result_list$beta_p[-wexc]*(result_list$beta_mu[-wexc]^2+result_list$beta_sigma[-wexc]))
		result_list2$z <- result_list$beta_chi[-wexc];
    result_list2$pval <- pchisq(result_list2$z^2,1,lower.tail=FALSE);
    result_list2$post <- result_list$beta_p[-wexc];
    result_list2$l0 <- result_list$l0_path;
    result_list2$modelEntropy <- -sum(modProb*log(modProb));
    result_list2$modelProb <- modProb;
    if(cleanSolution){
      fastlm <- function(y,X){
        n1 <- nrow(X)
        X <- cbind(rep(1,n1),X);
        ginv <- solve(t(X)%*%X);
        Xhat <- ginv%*%t(X);
        betahat <- Xhat%*%y;
        sig <- mean((y-X%*%betahat)^2)*((n1)/(n1-ncol(X)));
        zval <- betahat/sqrt(sig*diag(ginv));
        return(zval[-1]);
      }
      signif <- result_list2$pval < 0.05/m;
      #a1 <- fastlm(y,X[,-wexc][,signif])
      #print(a1)
      if(sum(signif)>1){
        result_list2$z[signif] <- fastlm(y,X[,-wexc][,signif]);
        result_list2$pval[signif] <- pchisq(result_list2$z[signif]^2,1,lower.tail=FALSE);
      }
    }
    return(result_list2);
	}else{
    result_list2 <- list();
    result_list2$beta <- result_list$e_beta[-wexc,];
    result_list2$alpha <- result_list$e_beta[wexc,];
    result_list2$z <- result_list$beta_chi[-wexc,];
    result_list2$pval <- pchisq(result_list2$z^2,1,lower.tail=FALSE);
    result_list2$post <- result_list$beta_p[-wexc,];
    result_list2$l0 <- result_list$l0_path;
    
    result_list2$modelEntropy <- -sum(modProb*log(modProb));
    result_list2$modelProb <- modProb;
    if(!is.null(result_list$kl_index)){
      result_list2$kl_index <- result_list$kl_index;
    }
    if(!is.null(result_list$kl)){
      result_list2$kl <- result_list$kl;
      result_list2$kl_min <- result_list$kl_min;
      result_list2$kl_mean <- result_list$kl_mean;
      result_list2$kl_se <- result_list$kl_se;
    }
        
    return(result_list2)
	}
}
