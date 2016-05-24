ppca.metabol <-
function(Y, minq=1, maxq=2, scale="none", epsilon = 0.1, plot.BIC=FALSE, printout=TRUE)
{

if (missing(Y)) {
        stop("Spectral data are required to fit the PPCCA model.\n")
}
if (minq > maxq) {
        stop("minq can not be greater than maxq.\n")
    }
if (maxq > ncol(Y)) {
        stop("maxq can not be greater than the number of variables.\n")
    }
if(maxq > 30)
{
	cat("Warning! Model fitting may become unstable for q > 30.\n\n")
}
if (epsilon > 1) {
        cat("Warning! Poor model covergence expected for epsilon > 1.\n")
    }
if (epsilon < 0.0001) {
        cat("Warning! Model covergence becomes very slow for epsilon < 0.0001.\n")
}


### Set up
V<-5000						 ## Maximum number of iterations
N<-nrow(Y)					 ## Number of observations
p<-ncol(Y)				     ## Number of variables
if(p > 375)
{
	stop("Spectral dimension is too large for current computation capabilities. Reduce the number of spectral bins to less than 375.")
}
				

### Storage for fitted models
Sig_q<-rep(0,maxq)                      ## Storage of noise covariance 
W_q<-list()                          ## Storage of loadings for Q models
U_q<-list()                          ## Storage of scores for Q models
AIC<-rep(0,maxq)                        ## Storage of AIC values
BIC<-rep(0,maxq)                        ## Storage of BIC values
ll<-matrix(0,maxq,V)                    ## Log likelihood for each iteration
lla<-matrix(0,maxq,V)                   ## Estimate of the asymptotic maximum of the log likelihood for each iteration
   
Y<-as.matrix(scaling(Y, type=scale))
##### Prior Parameters
Vp<-10                                               
C2p<-p*3
muhat<-colMeans(Y)					#  mu MLE
Yc<-sweep(Y,2,muhat,"-")                 ## Center data
S<-(1/nrow(Yc))*(t(Yc)%*%Yc)                   ## Empirical covariance matrix
temp<-eigen(S)   

### Fit maxq-PPCA Models 
for(q in minq:maxq)
{    
   ### Initialize the model 
   Sig<-abs((1/(p-q))*sum(temp$val[(q+1):p]))     ## Starting value for variance
   W<-temp$vec[,1:q]                              ## Starting value for loadings
   u<-t(cmdscale(dist(Y),q))

   tol <- epsilon+1                               ## Initial tolerance value
   v <- 0                                         ## Initializing iteration counter
   while(tol>epsilon)                             ## Until convergence
   {
      v <- v+1
      
      ## M step	  
      k<-S%*%W
      M_1<-solve((t(W)%*%W + Sig*diag(q)))  
      
      W<-k%*%solve(Sig*diag(q) + (M_1%*%t(W))%*%k)
     
      MLESig<-(1/p)*sum(diag(S - k%*%(M_1%*%t(W))))
      Sig<- c(((N*p)*MLESig + C2p)/((N*p) + Vp + 2))
            
      ## E step
      u<-M_1%*%(t(W)%*%t(Yc))

      ### Calculate observed log-likelihood 
      ll[q,v]<-sum(dmvnorm(Y, muhat, W%*%t(W)+Sig*diag(p), log=TRUE))
    
      ### Covergence assessment
      converge<-Aitken(ll, lla, v, q, epsilon)
      tol<-converge[[1]]
      lla[q,v]<-converge[[2]]
      
      if(v == V)
      { 
      	if(printout == TRUE)
   		{
      		cat("Algorithm stopped for q = ", q,". Maximum number of iterations exceeded.\n\n")
      	 }
      	tol<-epsilon-1
      }
   } # close while loop for PPCA

   if(printout == TRUE)
   {
   	cat("q = ", q, ": PPCA converged.\n\n")
   }

   ### Evaluate model selection criteria
   params<-(p*q) - (0.5*q*(q-1)) + 1		## Number of model parameters
   AIC[q]<-(2*ll[q,v]) - (2*params)
   BIC[q]<-(2*ll[q,v]) - (params*log(N))
      
   ####### Store model fit  
   U_q[[q]]<-u 
   Sig_q[q]<-Sig 
   W_q[[q]]<-W
}# Close q loop

#### Report optimal model
qopt<-c(minq:maxq)[BIC[minq:maxq]==max(BIC[minq:maxq])]
Uopt<-t(U_q[[qopt]])
Wopt<-W_q[[qopt]]
Sigopt<-Sig_q[qopt]

if(plot.BIC == TRUE)
{
	plot(minq:maxq, BIC[minq:maxq], type="b", xlab="q", ylab="BIC", col.lab="blue")
	abline(v=qopt,col="red", lty=2)
}

list(q=qopt, sig=Sigopt, scores=Uopt, loadings=Wopt, BIC=BIC[minq:maxq], AIC=AIC[minq:maxq])
}

