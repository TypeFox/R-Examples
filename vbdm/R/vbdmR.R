#R implementation of discrete mixture regression model (vanilla) for rare variant association analysis
#Benjamin A. Logsdon, Copyright 2012
#equations that are referred to below are in a companion discrete_mixture.pdf file

#############################discrete.mixture.test########################
#INPUTS:
#y: continous outcome to be tested
#G: A matrix of genotypes
#X: A matrix of confounding covariates
#tol: Optionally specified tolerance for convergence
#thres: Minor allele frequency threhold for variants within G
#scaling: Whether or not columns of G should be scaled to have mean 0 and variance 1
#hyper: Hyper parameters alpha and beta for prior over the mixing pprobability parameter.  Default is 2 and 2.
##########################################################################
#OUTPUTS:
#keep: the columns of G that passed the MAF threshold and are used in the analysis
#y: outcome
#gssq: the sum of squares of columns of G
#Gq: The unscaled G matrix
#X: the covariate matrix used in the analysis
#n: sample size
#m: number of features
#p: number of covariates
#pvec: posterior probability of association of each variant
#prob: the mixing probability (e.g. proporition of associated variants in a gene)
#theta: the estimate of the effect of the burden of associated variants
#gamma: the fixed effects
#sigma: the estimated residual variance
#Gp: G%*%pvec; ie the matrix vector product between genotypes and posterior probabilities
#resid: The residuals from the model
#lb: the lower bound of the alternative model
#lb_null: the lower bound of the null model
#tol: the tolerance used
#lrt: the likelihood ratio test
#p.value: the p-value of the association (assuming the lrt ~ chi^2_1)
##########################################################################
#EXAMPLE:
#set.seed(1)
#n <- 2000;
#m <- 20;
#p <- 2;
#G <- matrix(rbinom(n*m,2,0.01),n,m);
#y <- G%*%rbinom(m,1,.2);
#y <- y+ rnorm(n,0,sqrt(10));
#X <- matrix(rnorm(n*p),n,p);
#res_scaled <- discrete.mixture.test(y=y,G=G,X=X,scaling=T);
#res_unscaled <- discrete.mixture.test(y=y,G=G,X=X,scaling=F);

vbdmR <- function(y,G,X=NULL,tol=1e-4,thres=0.05,scaling=TRUE,hyper=c(2,2)){
  update_p <- function(model){
	  #This function produces the estimates of the posterior probabilities for
	  #each variant and stores them in the list element model$pvec.
  	#It also updates the overall residuals and probability
    #weighted burden scores after accounting for
    #each variant in the model when weighted by its posterior
    #probability.
    for(j in 1:model$m){
      #extract the current value of the posterior probability
      pold <- model$pvec[j];
      
      #compute the dot product between the residual *without* variant j
      #and the genotype vector for variant j
      vec1 <- t(model$resid)%*%model$G[,j]+model$gssq[j]*pold*model$theta;
      
      #compute the first term in the expression for the posterior probability
      #of variant j(related to the standard error of the estimate)
      #equation 10
      a1 <- as.numeric(model$theta^2*model$gssq[j]);
      #cat(a1,'a1\n')
      
      #compute the second term in the expression for the posterior probability
      #of variant j (related to the estimate of the effect of variant j)
      #equation 10
      a2 <- -as.numeric(2*model$theta*(vec1));
      #cat(a2,'a2\n')
      
      #compute the third term in the expression for the posterior probability
      #of variant j related to the overall mixture probability of being non-zero.
      #The higher this is the weaker the evidence of association needs to be to have a high
      #posterior probability
      #a3 <- -log(model$prob);
      #equation 10
      a3 <- -(digamma(model$prob*(model$m+sum(model$hyper)))-digamma(model$m+sum(model$hyper)));
      #cat(a3,'a3\n')
      #compute the fourth term in the expression for the posterior probability
      #of variant j (related to the overall mixture probability of being zero.
      #the higher this is, the stronger the evidence of association needs to be to have
      #a high posterior probability
      #a4 <- log(1-model$prob);
      #equation 10
      a4 <- (digamma((1-model$prob)*(model$m+sum(model$hyper)))-digamma(model$m+sum(model$hyper)))
      #cat(a4,'a4\n')
      #combine all the terms (essentially rescaling a1 and a2 by the relative scale of the
      #residuals so that all variants are treated on the same scale.
      #equation 10
      a5 <- (1/(2*model$sigma))*(a1+a2)+a3+a4;
      #cat(a5,'a5\n')
      
      #compute the logistic function of a5 to transform it from a real number 
      #into the posterior probability
      #equation 10
      
      model$pvec[j] <- 1/(1+exp(a5));
      #cat(model$pvec[j],'pj\n');
      
      #update the model residuals by adding the old effect and subtracting the new effect
      model$resid <- model$resid + model$theta*model$G[,j]*(pold-model$pvec[j]);
      
      #update the reweighted burden score, model$Gp, by adding the new probability reweighted
      #score, and subtracing the old probability reweighted score
      model$Gp <- model$Gp + model$G[,j]*(model$pvec[j]-pold);
      
    }
    
    #update the estimated of the overall probabilty of being non-zoer, model$prob:
    #equation 17/18
    model$prob <- (sum(model$pvec)+model$hyper[1])/(model$m+sum(model$hyper));
    #cat(model$prob,'prob\n');
    
    return(model);	
  }
  
  
  
  update_theta_gamma <- function(model){
    
    #store old value of theta
    theta_old <- model$theta;
    
    
    #compute numerator of new theta estimate
    #equation 12
    theta_new <- t(model$resid)%*%model$Gp + t(model$Gp)%*%model$Gp*model$theta;
    #cat(theta_new,'theta_new\n')
    #compute denominator of new theta estimate
    #equation 12
    const2 <- sum(model$Gp^2)+t(model$gssq)%*%(model$pvec-model$pvec^2);
    #cat(const2,sum(model$Gp^2),model$gssq,model$pvec,'denom\n')
    #update theta
    #equation 12
    theta_new <- as.numeric(theta_new/as.numeric(const2));
    #cat(theta_new,'tn2\n')
    if(!is.null(model$theta2)){
      #update theta under null model
      theta_new <- model$theta2;
    }
    model$theta<-theta_new;
    
    #update residuals
    model$resid <- model$resid + model$Gp*(theta_old-theta_new);
    
    #extract covariates
    Z <- model$X;
    #store old value of covariate effects
    #equation 14
    alpha_old <- c(model$gamma);
    
    
    #solve for new value of covariate effects
    #equation 14
    alpha <- solve(t(Z)%*%Z)%*%t(Z)%*%(model$y-model$Gp*theta_new);
    #cat(alpha,'alpha\n')
    
    #update the covariate effects
    #equation 14
    model$gamma <- alpha;
    
    #update the residuals
    model$resid <- model$resid + Z%*%(alpha_old-alpha);
    
    return(model);
  }
  update_sigma <- function(model){
    
    
    #compute correction for sum of square errors from lower bound
    #equation 15
    const2 <- model$gssq%*%(model$pvec-model$pvec^2);
    
    #compute error variance sigma including the correction
    #equation 15
    model$sigma <- (t(model$resid)%*%model$resid) + const2*model$theta^2;
    #cat(model$sigma,'sigma\n')
    
    #update sigma
    #equation 15
    model$sigma <- model$sigma/model$n;
    #cat(model$sigma,'sigma2\n')
    return(model);
    
    
    
  }
  update_lb <- function(model){
    
    
    #compute likelihood portion of lower bound
    #equation 16
    lb <- -0.5*(model$n*(log(2*pi*model$sigma)+1));
    ##cat(lb,"ll\n")
    #lba <- lb;
    #model$lb <- lb;
    
    alpha1 <- sum(model$pvec)+model$hyper[1];
    beta1  <- model$m-sum(model$pvec)+model$hyper[2];
    
    
    #the prior portion of the lower bound
    #equation 19
    lb <- lb + (digamma(alpha1)-digamma(alpha1+beta1))*(sum(model$pvec)+1);
    lb <- lb + (digamma(beta1)-digamma(alpha1+beta1))*(model$m-sum(model$pvec)+1);
    ##cat(lb,'elp\n')
    
    #entropy portion of lower bound
    #equation 20
    ent_1 <- model$pvec*log(model$pvec);
    ent_0 <- (1-model$pvec)*log(1-model$pvec);
    ent_1[is.na(ent_1)]<-0;
    ent_0[is.na(ent_0)]<-0;
    lb <- lb - sum(ent_1);
    lb <- lb - sum(ent_0);
    ##cat(model$pvec,'prob\n')
    ##cat(model$theta,'theta\n')
    ##cat(lb,'hz\n')
    lb <- lb + lbeta(alpha1,beta1) - (alpha1-1)*digamma(alpha1)
    lb <- lb - (beta1-1)*digamma(beta1) + (alpha1+beta1-2)*digamma(alpha1+beta1);
    ##cat(lb,'hp\n')
    
    model$lb <- lb;
    return(model);
  }
  
  run_discrete_mixture_vanilla <- function(y,G,X=NULL,tol=1e-4,theta=NULL,thres=0.05,scaling=TRUE,hyper=c(2,2)){
    model <- list();
    #preprocess covariates, thres is 
    if(is.null(X)){
      X <- as.matrix(rep(1,nrow(G)));
    }else{
      X <- cbind(rep(1,nrow(X)),X);
    }
    #preprocess G
    cs <- colMeans(G,na.rm=T)/2;
    cvec <- apply(cbind(cs,1-cs),1,min);
    keep <- which((cvec<thres)*(cvec>0)==1);
    cvec2 <- apply(cbind(cs,1-cs),1,which.min);
    G <- as.matrix(G[,keep]);
    if(length(keep)==0){
      stop("No rare variants left.\n");
    }
    flip <- which(cvec2[keep]==2);
    if(length(flip)>0){
      G[,flip] <- 2-G[,flip];
    }
    nav <- which(is.na(G));	
    if(length(nav)>0){
      G[nav]<-0;
    }
    Gq <- G;
    
    #scale all variants to have mean 0 and variance 1
    if(scaling){
      G <- as.matrix(scale(G));
    }else{
      G <- as.matrix(G);
    }
    #initialize model;
    model$keep <- keep;
    model$y <- y;
    model$G <- G;
    model$Gq <- Gq;
    model$gssq <- colSums(G^2);
    model$X <- X;
    model$n <- length(y);
    model$m <- ncol(G);
    model$p <- ncol(X);
    model$hyper <- hyper;
    
    model$pvec <- rep(0.5,model$m);
    model$prob <- 0.5;
    model$theta <- 0;
    model$gamma <- rep(0,model$p);
    model$sigma <- var(model$y);
    model$Gp <- model$G%*%model$pvec;
    model$resid <- model$y;
    model <- update_lb(model);
    model$tol <- tol;
    model$theta2 <- theta;
    lb_old <- -Inf;
    lb_new <- 0;
    
    while(abs(lb_old-lb_new)>model$tol){
      lb_old <- lb_new;
      model<-update_p(model);
      #print(model$pvec)
      model<-update_theta_gamma(model);
      model<-update_sigma(model);
      model<-update_lb(model);
      lb_new <- model$lb;
      #print(lb_new);
    }
    #print(model$pvec);
    #print(model$lb);
    return(model);
    
  }


	model_a <- run_discrete_mixture_vanilla(y=y,G=G,X=X,tol=tol,thres=thres,scaling=scaling,hyper=hyper);
	model_b <- run_discrete_mixture_vanilla(y=y,G=G,X=X,tol=tol,theta=0,thres=thres,scaling=scaling,hyper=hyper);
	model_a$lrt <- -2*model_b$lb + 2*model_a$lb;
	model_a$p.value <- pchisq(model_a$lrt,1,lower.tail=F);
	model_a$lb_null <- model_b$lb
	return(model_a);
}





