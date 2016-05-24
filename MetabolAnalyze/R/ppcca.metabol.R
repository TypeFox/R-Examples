ppcca.metabol <-
function(Y, Covars, minq=1, maxq=2, scale="none", epsilon = 0.1, plot.BIC=FALSE, printout=TRUE)
{
Y<-as.matrix(Y)
Covars<-as.matrix(Covars)
if (missing(Y)) {
        stop("Spectral data are required to fit the PPCCA model.\n")
    }
if (missing(Covars)) {
        stop("Covariate data are required to fit the PPCCA model.\n ")
    }
if (nrow(Y) != nrow(Covars)) {
        stop("Spectral data and covariate data should have the same number of rows.\n")
    }
if (missing(minq)) {
        minq<- 1
    }
if (missing(maxq)) {
        maxq<- 2
    }
if (minq > maxq) {
        stop("minq can not be greater than maxq.\n")
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


### Set up
V<-4000                       ## Maximum number of iterations
N<-nrow(Y)					  ## Number of observations
p<-ncol(Y)					  ## Number of variables
L<-ncol(Covars)                    ## Number of covariates
Covars<-standardize(Covars)			  ## Standardize covariates for stability
Covars<-rbind(rep(1, N), t(Covars))

### Checking the dimensionality of the data
if(p > 375)
{
	stop("Spectral dimension is too large for current computation capabilities. Reduce the number of spectral bins to less than 375.")
}


######## Storage Q-PPCA Models #########
Sig_q<-rep(0,maxq)                      ## Storage of noise covariance 
W_q<-list()                             ## Storage of loadings for Q models
Alpha_q<-list()                         ## Jackknife coefficients storage for Q models
U_q<-list()                             ## Storage of scores for Q models
AIC<-rep(0,maxq)                        ## Storage of AIC values
BIC<-rep(0,maxq)                        ## Storage of BIC values
ll<-matrix(NA,maxq,V)                    ## Log likelihood for each iteration
lla<-matrix(NA,maxq,V)                   ## Estimate of the asymptotic maximum of the log likelihood on each iteration
   
   
#unscaledY<-Y								
Y<-as.matrix(scaling(Y, type=scale))           ## Scale the data
##### Prior Parameters
Vp<-10                                               
C2p<-p*3
muhat<-colMeans(Y)					#  mu MLE
Yc<-sweep(Y,2,muhat,"-")                 ## Center data
S<-(1/nrow(Yc))*(t(Yc)%*%Yc)                   ## Empirical covariance matrix
temp<-eigen(S)    
   
### Fit maxq-PPCCA Models 
for(q in minq:maxq)
{    

   ### Initializing the model 
   Sig<-abs((1/(p-q))*sum(temp$val[(q+1):p]))   ## Starting values for variance
   W<-temp$vec[,1:q]                            ## Starting value for loadings
   scores<-t(solve((t(W)%*%W) + (Sig*diag(q)))%*%t(W)%*%t(Yc))
   Alpha<-matrix(0, q, L+1)
   for(i in 1:q)
   {
   	if(L==1)
   	{
   		dat<-data.frame(cbind(scores[,i], as.matrix(Covars[2:(L+1),])))
      }else{
            dat<-data.frame(cbind(scores[,i], t(Covars[2:(L+1),])))
      }
   	Alpha[i,]<-glm(dat, family=gaussian)$coef
   }
      
   tol <- epsilon+1                                          ## Initial tolerance value
   v <- 0                                                    ## Initializing iteration counter
   while(tol>epsilon)                                        ## Until convergence
   {
      v <- v+1
      
      # E-step      
      M_1<-solve(t(W)%*%W + Sig*diag(q))                         
      u<-M_1%*%(t(W)%*%t(Yc) + Sig*(Alpha%*%Covars))   ## Expected value of scores.
      Sum_Euu<-(nrow(Yc)*Sig*M_1) + (u%*%t(u))
      
      ## M step
      Alpha<-(u%*%t(Covars))%*%solve(Covars%*%t(Covars))        ## Estimation of regression coefficients
      
      W<-(t(Yc)%*%t(u))%*%solve(Sum_Euu)         ## Estimation of loadings.                         
      
      YWEu<-sum(diag(Yc%*%W%*%u))
      MLESig<-(nrow(Yc)*sum(diag(S)) + sum(diag((t(W)%*%W)%*%Sum_Euu)) - 2*YWEu)/(p*nrow(Yc))    ## Estimation of the variance MLE.
      Sig<- c(((N*p)*MLESig + C2p)/((N*p) + Vp + 2))

      ### Calculate observed log-likelihood 
      Den<-rep(NA, nrow(Y))
      Sigma<-W%*%t(W)+(Sig*diag(p))
      mumat<-W%*%(Alpha%*%Covars) + matrix(muhat, nrow=p, ncol=N, byrow=FALSE)    
      for(i in 1:nrow(Y))
      {
         Den[i]<-(dmvnorm(Y[i,], mumat[,i], Sigma, log=TRUE))
      }
      ll[q,v]<-sum(Den)

      ### Covergence Assessment 
      converge<-Aitken(ll, lla, v, q, epsilon)
      tol<-converge[[1]]
      lla[q,v]<-converge[[2]]
      if(v == V)
      { 
      		cat("Algorithm stopped for q = ", q,". Maximum number of iterations exceeded.\n\n")
      		tol<-epsilon-1
      }
      } #while loop for PPCCA
     
   if(printout == TRUE)
   {
   	cat("q = ", q, ": PPCCA converged.\n\n")
   }

   ####### Evaluate model selection criteria
   params<-(p*q) - (0.5*q*(q-1)) + (q*(L+1)) + 1		## Number of model parameters
   AIC[q]<-2*ll[q,v] - (2*params)
   BIC[q]<-2*ll[q,v] - params*log(N)
      
   ####### Store model fit  
   U_q[[q]]<-u 
   Sig_q[q]<-Sig 
   W_q[[q]]<-W
   Alpha_q[[q]]<-Alpha
   }#q

#### Report optimal model
qopt<-c(minq:maxq)[BIC[minq:maxq]==max(BIC[minq:maxq])]
Uopt<-t(U_q[[qopt]])
Wopt<-W_q[[qopt]]
Sigopt<-Sig_q[qopt]
Alphaopt<-Alpha_q[[qopt]]

if(plot.BIC == TRUE)
{
	plot(minq:maxq, BIC[minq:maxq], type="b", xlab="q", ylab="BIC", col.lab="blue")
	abline(v=qopt,col="red", lty=2)
}

list(q=qopt, sig=Sigopt, scores=Uopt, loadings=Wopt, coefficients=Alphaopt, BIC=BIC[minq:maxq], AIC=AIC[minq:maxq])
} # End ppcca.metabol

