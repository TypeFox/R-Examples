mppca.metabol <-
function(Y, minq=1, maxq=2, ming, maxg, scale="none", epsilon = 0.1, plot.BIC = FALSE)
{
if (missing(Y)) {
        stop("Spectral data are required to fit the PPCCA model.\n")
    }
if (missing(ming)) {
        stop("The user needs to specify the minimum number of groups to be fitted.\n")
    }
if (missing(maxg)) {
        stop("The user needs to specify the maximum number of groups to be fitted.\n")
    }
if (minq > maxq) {
        stop("minq can not be greater than maxq.\n")
    }
if (ming > maxg) {
        stop("ming can not be greater than maxg.\n")
    }
if (maxq > ncol(Y)) {
        stop("maxq can not be greater than the number of variables.\n")
    }
if(maxq > 10)
{
	cat("Warning! Model fitting may become very slow for q > 10.\n\n")
}
if (epsilon > 1) {
        cat("Warning! Poor model covergence expected for epsilon > 1.\n")
    }
if (epsilon < 0.0001) {
        cat("Warning! Model covergence becomes very slow for epsilon < 0.0001.\n")
    }
    

### Set up of values and storage
V<-40000                      ## Maximum number of iterations
N<-nrow(Y)					  ## Number of observations
p<-ncol(Y)					  ## Number of variables
reset<-FALSE                        

### Checking the dimensionality of the data
if(p > 375)
{
	stop("Spectral dimension is too large for current computation capabilities. Reduce the number of spectral bins to less than 375.")
}
if(round(N/maxg) < 5)
{
	cat("Warning! Model fitting may become unstable as maxg is large relative to N.\n\n")
}

Y<-as.matrix(scaling(Y, type=scale))           ## Scale the data


### STORAGE maxg*maxq-PPCA MODELS
W_g<-list()             ## Storage of all the loadings for all the models                     
Mu_g<-list()            ## Mean vectors for all the models
Pi_g<-list()            ## Prior probabilities of membership for all the models
Tau_g<-list()           ## Prior probabilities of membership for all the models
Sig_g<-list()           ## Covariance of the noise process for all the models
ll<-matrix(0,maxq,V)    ## Log likelihood for each iteration
lla<-matrix(0,maxq,V)   ## Estimate of the asymptotic maximum of the log likelihood on each iteration
AIC<-matrix(-Inf,maxg,maxq)
BIC<-matrix(-Inf,maxg,maxq)


### Fit maxg*maxq-PPCA MODELS 
for(g in ming:maxg)     
{
   Sig_q<-list()            ## Covariance of the noise for Q models 
   W_q<-list()              ## Loadings for Q models 
   Pi_q<-list()             ## Group membership prob for Q models 
   mu_q<-array(NA, c(p,g,maxq))  ## Group means for Q models
   Tau_q<-array(NA, c(nrow(Y),g,maxq)) ## Posterior probs for Q models
 
   ### FIT maxq-MPPCA MODELS WITH g GROUPS 
   for(q in minq:maxq)
   {                                                
         ### INITIALIZING THE MODEL                                                                         
         if(g == 1)
         { 
            Tau<-matrix(rep(1,nrow(Y))) 
            logTau<-Tau             
            Pi<-1
            mu<-matrix(colMeans(Y), p, 1)
         }else{                                 
			temp<-mclustBIC(cmdscale(dist(Y)), G=g)
			res<-summary(temp, cmdscale(dist(Y)))
			Tau<-res$z
            logTau<-Tau
            Pi<-res$parameters$pro  
            if(any(table(map(Tau)) < (q+1))|ncol(unmap(map(Tau)))<g)
	        {    
               Tau<-unmap(sample(rep(1:g, 2), nrow(Y), replace=TRUE))
               Pi<-apply(Tau,2,sum)/nrow(Y)
               logTau<-Tau
	        }#if     
            mu<-matrix(0,p,g)  
			for(i in 1:g)
			{   
		   		mu[,i]<-apply(Y[map(Tau)==i,], 2, mean)
			}      
          }#else 
                                                         
         W<-array(NA, c(p,q,g))                ## Loadings for MPPCA model
         sig<-rep(NA,g)                        ## Covariance of the noise process 
         for(i in 1:g)
         {
           S<-cov(Y[map(Tau)==i,])
           decomp<-eigen(S)
           W[,,i]<-decomp$vec[,1:q]
           sig[i]<-abs(mean(decomp$val[(q+1):p]))
           if(sig[i]<0.0001){sig[i]<-abs(rnorm(1,0.01, 0.001))}
         }#i
         Sig<-sum(sig)/g

         ### CONVERGENCE CRITERIA 
         tol <- epsilon+1    ## Initial tolerance value
         v <- 0              ## Initializing iteration counter

         while(tol>epsilon)  ## Until convergence
         {
            v<-v+1
            
            ## First E-step
            res<-estep1(Y, Tau, Pi, mu, W, Sig, g, p, reset)                            ## E-Step of the first cycle
            Tau<-res[[1]]
            logTau<-res[[2]]
            reset<-res[[3]] 

            ## First M-step
	        res1<-mstep1(Y, Tau, Pi, mu, g)                                                       ## M-Step of the first cycle
            Pi<-res1[[1]] 
            mu<-res1[[2]]

            ## Second E-step
            res<-estep2(Y, Tau, Pi, mu, W, Sig, g, p, reset)                                        ## E-Step of the second cycle 
            Tau<-res[[1]]
            logTau<-res[[2]]
            reset<-res[[3]]

            ## Second M-step
	        res2<-mstep2(Y, Tau, Pi, mu, W, Sig, g, p, q)                                     
    	    W<-res2[[1]]
	        Sig<-res2[[2]]

            ### OBSERVED LOG LIKELIHOOD
            store<-matrix(NA,nrow(Y),g)
            for(i in 1:g)
            {
               store[,i]<-Pi[i]*dmvnorm(Y, mu[,i], W[,,i]%*%t(W[,,i])+Sig*diag(p))
             }#g
            if(reset == FALSE)
            {
            	ll[q,v]<-sum(log(rowSums(store)))
            }else{
	            ll[q,v]<-sum(apply(logTau, 1, max) + log(apply(exp(logTau-apply(logTau,1,max)) ,1, sum)))
	            reset<-FALSE
	      }


            ### CONVERGENCE ASSESSMENT 
            Converg<-Aitken(ll, lla, v, q, epsilon)
            tol<-Converg[[1]]
            lla[q,v]<-Converg[[2]]

            if(v == V)
            { 
      	    cat("Algorithm stopped for g = ", g, "q =", q, ". Maximum number of iterations exceeded.\n\n")
      	    tol<-epsilon-1
            }
         }#while loop
         cat("g = ", g, ",", "q = ", q, ": MPPCA converged.\n\n")

         ### EVALUATING BEST FIT MODEL
         params<-(g-1) + (g*p) + (g*((p*q) - ((q*(q-1))/2))) + 1                         ## Number of model parameters
         AIC[g,q]<-(2*ll[q,v]) - (2*params) 
         BIC[g,q]<-(2*ll[q,v]) - (params*log(nrow(Y)))

         ### STORAGE Q-PPCA MODELS 
         Tau_q[,,q]<-Tau
         Sig_q[[q]]<-Sig
         mu_q[,,q]<-mu
         W_q[[q]]<-W
         Pi_q[[q]]<-Pi

         ### STORAGE G*Q-PPCA MODELS 
         Sig_g[[g]]<-Sig_q
         Tau_g[[g]]<-Tau_q
         W_g[[g]]<-W_q
         Pi_g[[g]]<-Pi_q
         Mu_g[[g]]<-mu_q

     }#q loop
}#g loop

### BEST MODEL 
qgopt<-which(BIC==max(BIC), arr.ind=TRUE)
gopt<-qgopt[1,1]
qopt<-qgopt[1,2]

Wopt<-W_g[[gopt]][[qopt]]
Piopt<-Pi_g[[gopt]][[qopt]]
Muopt<-Mu_g[[gopt]][,,qopt]
Sigopt<-Sig_g[[gopt]][[qopt]]
Tauopt<-Tau_g[[gopt]][,,qopt]


######### SCORES FOR EACH GROUP 
Uopt<-list()
for(g in 1:gopt)
{
  M_1<-solve((t(as.matrix(Wopt[,,g]))%*%as.matrix(Wopt[,,g])) + ((Sigopt)*diag(qopt)))
  if(gopt==1)
  {
  	Uopt[[g]]<-t((M_1)%*%(t(as.matrix(Wopt[,,g]))%*%t(sweep(Y, 2, Muopt, "-"))))
  }#if
  else{
  	Uopt[[g]]<-t((M_1)%*%(t(as.matrix(Wopt[,,g]))%*%t(sweep(Y[map(Tauopt)==g,], 2, Muopt[,g], "-"))))
  	}#else
  	colnames(Uopt[[g]])<-paste("PC_",1:qopt, space="")
}#end

if(minq==maxq|ming==maxg)
{
if(minq==maxq & ming!=maxg)
{
	BIC<-matrix(BIC[ming:maxg, minq:maxq],maxg-ming+1, maxq-minq+1, byrow=TRUE)
	AIC<-matrix(AIC[ming:maxg, minq:maxq],maxg-ming+1, maxq-minq+1, byrow=TRUE)
}
if(ming==maxg & minq!=maxq)
{
	BIC<-matrix(BIC[ming:maxg, minq:maxq],maxg-ming+1, maxq-minq+1, byrow=FALSE)
	AIC<-matrix(AIC[ming:maxg, minq:maxq],maxg-ming+1, maxq-minq+1, byrow=FALSE)
}
if(ming==maxg & minq==maxq)
{
	BIC<-as.matrix(BIC[ming:maxg, minq:maxq])
	AIC<-as.matrix(AIC[ming:maxg, minq:maxq])
}
}else{
	BIC<-as.matrix(BIC[ming:maxg, minq:maxq])
    AIC<-as.matrix(AIC[ming:maxg, minq:maxq])
}
rownames(BIC)<-paste("G =", ming:maxg)
colnames(BIC)<-paste("q =", minq:maxq)
rownames(AIC)<-rownames(BIC)
colnames(AIC)<-colnames(BIC)


### BIC plot
if(plot.BIC == TRUE)
{
   if((maxg-ming) == 0)
   {
      plot(minq:maxq, BIC, type="b", xlab="q", ylab="BIC", col.lab="blue")
      abline(v=qopt,col="red", lty=2)
   }
   if((maxq-minq) == 0)
   {
      plot(ming:maxg, BIC, type="b", xlab="G", ylab="BIC", col.lab="blue")
      abline(v=gopt,col="red", lty=2)
   }
   if( (ming!=maxg) & (minq!=maxq))
   {
      ht(matrix(BIC, maxg-ming+1, maxq-minq+1, byrow=FALSE), q=qopt, g=gopt, xlab = "PCs", ylab = "Groups",main = "BIC Values",Rowv = NA, Colv = NA, scale="none")
   }
}


if(gopt==1)
{
   clustering=Tauopt
}else{
   clustering=map(Tauopt)
}

list(q=qopt, g=gopt, sig=Sigopt, scores=Uopt, loadings=Wopt, Pi=Piopt, mean=Muopt, tau=Tauopt, clustering=clustering, BIC=BIC, AIC=AIC)
} # End mppca.metabol function

