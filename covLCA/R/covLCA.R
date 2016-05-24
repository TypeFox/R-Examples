covLCA <-
function(formula1,formula2,data,nclass=2,maxiter=1000,tol=1e-10,
                beta.start=NULL,alpha.start=NULL,gamma.start=NULL,beta.auto=TRUE,alpha.auto=TRUE,gamma.auto=TRUE,nrep=1,verbose=TRUE,calc.se=TRUE) #A: formula1 for LC probabilities, formula2 for conditional probabilities. The left-hand side of both formulae must be the same
				#A: values provided by the user must be correct; no check

{ 
	starttime <- Sys.time() #A : current date and time
	mf1 <- model.response(model.frame(formula1,data,na.action=NULL)) #A: manifest variables in formula 1
	mf2 <- model.response(model.frame(formula2,data,na.action=NULL)) #A: manifest variables in formula 2
	 
	if (any((as.numeric(mf1)-as.numeric(mf2))!=0,na.rm=TRUE))#A: check whether both manifest tables are the same
	{
		stop("\n ALERT: manifest variables in both formulae must be identical. \n \n")
		ret <- NULL
	}

	#A: given that mf1 and mf2 are identical, check whether some manifest variables values are either <1 or non-integer.
    if (any(mf1<1,na.rm=TRUE) | any(round(mf1) != mf1,na.rm=TRUE)){ 

	cat("\n ALERT: some manifest variables contain values that are not 
  positive integers.  For covLCA to run, please recode categorical  
  outcome variables to increment from 1 to the maximum number of 
  outcome categories for each variable. \n \n")
        ret <- NULL}else #A : if all values for manifest variables have the right format
	{ 
		mframe1 <- model.frame(formula1,data,na.action=NULL) #A: dataframe with all necessary columns for formula1. 
		miss1=is.na(mframe1) #A: missing values for latent class covariates
		ind.miss1=(apply(miss1,1,sum)>0) #A: cases with missing values for latent class covariates 
		mframe2 <- model.frame(formula2,data,na.action=NULL) #A: dataframe with all necessary columns for formula2.
		miss2=is.na(mframe2) #A: missing values for conditional covariates
		ind.miss2=(apply(miss2,1,sum)>0) #A: cases with missing values for conditional covariates
		
		mframe1=mframe1[!(ind.miss1)& !(ind.miss2),] #A: in covariates x, we remove cases with missing values either in x or in z
		mframe2=mframe2[!(ind.miss1)& !(ind.miss2),] #A: in covariates z, we remove cases with missing values either in x or in z
		
		y <- model.response(mframe1) #A: manifest variables
		if (any(sapply(lapply(as.data.frame(y),table),length)==1)) #A: 1) create a contingency table for each manifest variables, 2) check whether each contingency table contains only one cell
		{
			y <- y[,!(sapply(apply(y,2,table),length)==1)] #A : 1) create a contingency table for each manifest variable, 2) remove, from the manifest variables database, columns (variables) taking on only one value. 
			cat("\n ALERT: at least one manifest variable contained only one outcome category, and has been removed from the analysis. \n \n")
		}
		x <- model.matrix(formula1,mframe1) #A: matrix with 1st column of 1s, and then one column for each latent class covariate
		z <- model.matrix(formula2,mframe2) #A: matrix with 1st column of 1s, and then one column for each conditional covariate
		if (ncol(z)==2){z <- array(z[,2],dim=c(dim(z)[1],1))}else {z <- z[,2:dim(z)[2]]} #A: Don't need an intercept because gammas are intercepts
		
		N <- nrow(y) #A: number of individuals
		J <- ncol(y) #A: number of manifest variables
		K.j <- t(matrix(apply(y,2,max))) #A: vector containing, for each manifest variable, the number of categories.
		
		#A: We add a restriction to simplify: all manifest variables have the same number of categories
		if (length(unique(as.vector(K.j)))>1)
		{
			cat("\n ALERT: all manifest variables must have the same number of outcome categories. \n \n")
			ret <- NULL
		}
		#A: end of restriction
		
		R <- nclass #A: number of LC
		S1 <- ncol(x) #A: number of covariates for LC probabilities (including the intercept)
		S2=ncol(z) #A number of covariates for conditional probabilities (without the intercept)
		
	
		
		eflag <- FALSE
		probs.start.ok <- TRUE
		ret <- list() #A: list which will contain all results
		
		ret$llik <- -Inf
		ret$attempts <- NULL
		
		for (repl in 1:nrep) 
		{ # automatically reestimate the model multiple times to locate the global max llik
			error <- TRUE; firstrun <- TRUE #A: put to FALSE at the end of the "while(error)" loop
			
			bet <- beta.init <- beta.start
			alph <- alpha.init <- alpha.start
			gamm <- gamma.init <- gamma.start
			if (beta.auto)
			{
				bet <- covLCA.initialBeta(y,R,x) #A: Output : a vector of initial values for betas, beta_jp, jp=11,...,1P, 21,...,2P, (J-1)1,...,(J-1)P
				beta.initAuto <- bet
				
			}
			if (alpha.auto)
			{
				alphgamm=covLCA.initialAlphaGamma(y,z,R,K.j,S2)
				#A: output for alpha: a matrix where each row correspond to a manifest variable, and row m contains alpha_mpk, where pk=11,12,...,1(K-1),21,22,...,2(K-1),...,P1,P2,...P(K-1)
				#A: output for gamma: a matrix where each row correspond to a manifest variable, and row m contains a repetition of gamma_mk: m1,m2,...,m(K-1),m1,m2,...,m(K-1),...
				alph <-  alphgamm$Alpha #A: row m: alpha_mpk
				alpha.initAuto <- alph
			}
			if (gamma.auto)
			{
				gamm <-  alphgamm$Gamma #A: row m: gamma_m(j)k
				gamma.initAuto <- gamm
			}
			
			
			while (error) 
			{ # error trap
				error <- FALSE
				
				#A: if no initial parameters, or if re-start the while(error) loop, or if not 1st replication: generate initial parameters
				if ((is.null(beta.start)&!beta.auto) | (!firstrun) | (repl>1))
				{
					bet <- beta.init <- rnorm(S1*(R-1)) #A: random vector of length (nb covariates LC*(nb LC-1)) => input for covLCA.updatePrior
					
				}
				if ((is.null(alpha.start)& !alpha.auto) | (!firstrun) | (repl>1))
				{ 
					alph <- alpha.init <- array(data=rnorm(J*S2*(K.j[1]-1)),dim=c(J,S2*(K.j[1]-1))) #Adapted to "S2=1"
				}
				if ((is.null(gamma.start)& !gamma.auto) | (!firstrun) | (repl>1))
				{ 
					gamm <- gamma.init <- array(data=rnorm(J*K.j[1]*R), dim=c(J,(K.j[1]-1)*R)) #Adapted to "S2=1"
				} #A: otherwise : has assigned values given by the user
				
				
							
				prior <- covLCA.updatePrior(bet,x,R)  #A: compute LC probabilities based on the last value of the parameters
				probs <- covLCA.updateCond(alph,gamm,z,R,J,K.j,S2,N) #A: compute conditional probabilities based on the last value of the parameters 
				
				
								
				iter <- 1
				llik <- matrix(NA,nrow=maxiter,ncol=1) #A: column vector with as many elements as maxiter
				llik[iter] <- -Inf #A: first element of llik is -inf
				dll <- Inf #A: difference between the log-lik of 2 iterations
				
				
				
				while ((iter <= maxiter) & (dll > tol) & (!error)) #A: 2nd loop, to repeat if has not reached max number of iterations and the improvement is high enough and there is no error
				{
					iter <- iter+1
					cat("Iteration number ",iter,"\n")
					flush.console()
					#print(iter)
					cat("\n")
					
				#A: E step:	
					rgivy <- covLCA.postClass(prior,probs,y,K.j) #A: compute posterior probs ("R given y"), a matrix where rows=indiv, cols=LC
				
					
					
				#A: M step:	given the estimate for the posterior, look for alphas, betas and gammas which maximize function Q
											
					dd.bet <- covLCA.dQdBeta(rgivy,prior,x) #A: Gradient and Hessian of Q_beta wrt betas
					bet <- bet + ginv(-dd.bet$hess) %*% dd.bet$grad #A: update betas through the one-iteration Newton-Raphson method 
					prior <- covLCA.updatePrior(bet,x,R) #A: update LC probabilities 
			
					
					
					#A: update omegas (gammas and alphas)
					for (m in 1:J) #A: for each manifest variable (i.e., for each subset of parameters alphas and gammas)
					{
						dd.gam <- covLCA.dQdGamma(rgivy,probs,y,K.j,m) #A: Gradient and Hessian of Q_omega wrt gammas
						dd.alph <- covLCA.dQdAlpha(rgivy,probs,z,K.j,m,y,S2) #A: Gradient and Hessian of Q_omega wrt alphas
						dd.alph.gam <- covLCA.dQdAlphaGamma(rgivy,probs,z,K.j,m,S2) #A: Hessian of Q_omega wrt alphas and gammas
						hess1=cbind(dd.gam$hess,dd.alph.gam)
						hess2=cbind(t(dd.alph.gam),dd.alph$hess)
						hess=rbind(hess1,hess2) #A: Hessian of Q_omega wrt the full set of parameters omegas (alphas and gammas)
						
						new <- c(gamm[m,],alph[m,]) + ginv(-hess) %*% c(dd.gam$grad,dd.alph$grad) #A: update gammas and alphas
						gamm[m,] <- new[1:dim(gamm)[2]]
						alph[m,] <- new[(dim(gamm)[2]+1):length(new)]
					}
				
					#A: update conditional probabilities
					probs <- covLCA.updateCond(alph,gamm,z,R,J,K.j,S2,N)
					#cat("Conditional probabilities updated \n")
					
					#A: compute log-likelihood						
					llik[iter] <- sum(log(rowSums(prior*covLCA.ylik(probs,y,K.j)))) #A: sum_{i=1}^N * log(sum_{j=1}^J p_j*prod_{m=1}^M prod_{k=1}^K p_mjk^{y_imk}) = log-likelihood as in formula (13) in user guide of poLCA
					cat("Llik=",llik[iter],"\n")
					
					dll <- llik[iter]-llik[iter-1] #A: improvement in the log-likelihood
					cat("dll=",dll,"\n")
					if (is.na(dll))  #A: if the log-likelihood of this iteration is missing ? 
					{
						error <- TRUE
					}else if ((S1>1) & (dll < -1e-7)) #A: if has covariates and the improvement is actually a deterioration
					{
						error <- TRUE
					}
				}#A: end of 2nd loop
				
				
				rgivy2 <- covLCA.postClass(prior,probs,y,K.j) #A: compute posterior probabilities with the final parameter estimates
				
				if (!error) #A: if error=FALSE, ie significant improvement
				{ 
					if (calc.se) 
					{
						ParamVar <- covLCA.paramVariance(prior,probs,rgivy2,R,S1,S2,J,K.j,x,y,z,N) #A: matrix of Var-Cov of the parameters
					} 
					else 
					{	
						ParamVar <- NA
						#se <- list(probs=NA,P=NA,b=NA,var.b=NA)
					}
				} else #A: if has error
				if (error)
				{
					eflag <- TRUE
				}
				firstrun <- FALSE
			}
			# finish estimating model without triggering error
			#A: end of the loop "while(error)"
			
			
			ret$attempts <- c(ret$attempts,llik[iter])
			if (llik[iter] > ret$llik) #A: if last llik obtained is higher than the previous maximum: choose these estimates
			{
				ret$llik <- llik[iter]             # maximum value of the log-likelihood 
				
				ret$beta.start <- beta.init
				ret$alpha.start <- alpha.init
				ret$gamma.start <- gamma.init
				ret$beta.auto <- beta.auto       # A: logical, TRUE if user asked for automatic search of reasonable initial parameter values
				ret$alpha.auto <- alpha.auto
				ret$gamma.auto <- gamma.auto
				if(beta.auto) ret$beta.initAuto <- beta.initAuto
				if(alpha.auto)ret$alpha.initAuto <- alpha.initAuto
				if(gamma.auto)ret$gamma.initAuto <- gamma.initAuto
				
				ret$probs <- probs                     # A: Many many values !
				ret$prior <- prior
				
				ret$posterior <- rgivy             # NxR matrix of posterior class membership probabilities 
				ret$posterior2 <- rgivy2			# not in poLCA
				ret$predclass <- apply(ret$posterior,1,which.max)   # Nx1 vector of predicted class memberships, by modal assignment 
				ret$P <- colMeans(ret$posterior)   # estimated class population shares 
				names(ret$P) <- paste("Latent class",1:R,sep=" ")
				ret$numiter <- iter-1              # number of iterations until reaching convergence 
				
				if (S1>1) #A: if has covariates for LC probabilities
				{
					b <- matrix(bet,nrow=S1)  # A: matrix where rows=covariates, columns=LC
					rownames(b) <- colnames(x)
					colnames(b) <- paste(1:(R-1),"vs",R,sep=" ")
					
					ret$coeffBeta <- b                 # coefficient estimates (when estimated)
					ret$param.se <- sqrt(diag(ParamVar))
					ret$param.V <- ParamVar
					
				} 
				if (S2>=1) #A: if has covariates for conditional probabilities #Adapted to "S2=1"
				{
					g <- gamm
					rownames(g) <- colnames(y)  #A: put colnames !
					colnames(g) <- paste("LC ",rep(seq(1,R),rep((K.j[1]-1),R)),", k=",rep(1:(K.j[1]-1),R),sep="")
					ret$coeffGamma <- g
					a <- alph
					rownames(a) <- colnames(y)
					colnames(a) <- paste("Var. ",rep(seq(1,S2),rep((K.j[1]-1),S2)),", k=",rep(1:(K.j[1]-1),S2),sep="")
					ret$coeffAlpha <- a
					ret$meanProbs=covLCA.meanCond(a,g,z,J,K.j,R,S2,N) #A: conditional probabilities evaluated at the sample means of the covariates
					dimnames(ret$meanProbs)[[1]]=colnames(y)
					dimnames(ret$meanProbs)[[2]]=paste("Pr(",1:K.j[1],")",sep="")
					dimnames(ret$meanProbs)[[3]]=paste("Latent class",1:R,sep=" ")
				}
				
				ret$eflag <- eflag                 # error flag, true if estimation algorithm ever needed to restart with new initial values
			}
			if (nrep>1 & verbose) { cat("Model ",repl,": llik = ",llik[iter]," ... best llik = ",ret$llik,"\n",sep=""); flush.console() }
		} # end replication loop
		

		
		ret$npar <- S1*(R-1)+S2*J*(K.j[1]-1)+J*(K.j[1]-1)*R # number of degrees of freedom used by the model (number of estimated parameters)
		
		ret$aic <- (-2 * ret$llik) + (2 * ret$npar)         # Akaike Information Criterion 
		ret$bic <- (-2 * ret$llik) + (log(N) * ret$npar)    # Schwarz-Bayesian Information Criterion
		ret$Nobs <- sum(rowSums(y==0)==0)                   # number of fully observed cases (if na.rm=F)
		ret$identifiability <- covLCA.identifiability(J,K.j,R,x,z,ret$npar,ret$coeffAlpha,ret$coeffBeta,ret$coeffGamma)
	
		

		ret$y <- data.frame(y)             # outcome variables 
		ret$x <- data.frame(x)             # A: covariates for latent classes 
		ret$z <- data.frame(z)             # A: covariates for conditional probabilities 
		ret$N <- N                         # number of observations 
		ret$maxiter <- maxiter             # maximum number of iterations specified by user 
		if (ret$numiter==ret$maxiter) cat("ALERT: iterations finished, MAXIMUM LIKELIHOOD NOT FOUND \n \n") #A: from "print.poLCA.R"
		ret$resid.df <- min(ret$N,(prod(K.j)-1))-ret$npar # number of residual degrees of freedom 
		class(ret) <- "covLCA"                

		ret$time <- Sys.time()-starttime   # how long it took to run the model 
    
	}
	   return(ret)
}
