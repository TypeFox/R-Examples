ppcca.metabol.jack <-
function(Y, Covars, minq=1, maxq=2, scale="none", epsilon = 0.1, conflevel=0.95)
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
      Vj<-40000                    ## Maximum number of jackknife iterations
      N<-nrow(Y)			       ## Number of observations
      p<-ncol(Y)		        	  ## Number of variables
      Covars<-as.matrix(Covars, nrow=N, byrow=TRUE)
      L<-ncol(Covars)                      ## Number of covariates

	  ## Fit a PPCCA model to the full dataset
	  modelfit<-ppcca.metabol(Y, Covars, minq, maxq, scale=scale, epsilon=epsilon, plot.BIC=FALSE, printout=FALSE)               
      q<-modelfit$q
	  Wj<-Wopt<-modelfit$loadings
	  Sigj<-Sigopt<-modelfit$sig
      Alphaj<-Alphaopt<-modelfit$coefficients
      Uopt<-modelfit$scores
      BIC<-modelfit$BIC
      AIC<-modelfit$AIC      

      Covars<-standardize(Covars)			 ## Standardize covariates for stability
      Covars<-rbind(rep(1, N), t(Covars))
      if(sum(rownames(Covars)== NULL) == 0)
      {
      	rownames(Covars)<-c(1:(L+1))
      	rownames(Covars)[1]<-"Intercept"
      	for(i in 2:(L+1))
      	{
        	rownames(Covars)[i]<-paste("Covariate", i, sep="")
        } #i
      }
      rownames(Covars)[1]<-"Intercept"
      
	cat("PPCCA fitted to all ", N, " observations. Optimal model has",q,"factors.\n\n")
	
	### Storage for jackknife output
	Sigstore<-rep(NA, N)
	Wstore<-array(0,c(p,q,N))
    Alphastore<-array(0,c(q,L+1,nrow(Y)))
	ll<-matrix(0,q,Vj)                    ## Log likelihood for each iteration
	lla<-matrix(0,q,Vj)                   ## Estimate of the asymptotic maximum of the log likelihood on each iteration
      
    ### Jackknife resampling 
    cat("Using the jackknife to estimate parameter uncertainty...\n \n")
	for(n in 1:N)
	{
  	  J<-Y[-n,] 			   	  ## Jackknife data set
      Covarsj<-Covars[,-n] 				  ## Jackknife covariate data set
	  J<-scaling(J, scale)                ## Scaling
	  Vp<-10                                               
	  C2p<-p*3
	  muhat<-colMeans(J)	
	  Jc<-sweep(J,2,muhat,"-")  ## Center jackknife data
  	  Sj<-(1/nrow(Jc))*(t(Jc)%*%Jc)

	  tol<-epsilon+1
      v<-0 
      while(tol>epsilon)
      {
     		v <- v+1  
            if(v == 1)
     		{
     			Wj<-Wopt
     			Sigj<-Sigopt
     			Alphaj<-Alphaopt
     			uj<-t(Uopt[-n,])
     		}
     		
     		# E-step 
            M_1<-solve(t(Wj)%*%Wj + Sigj*diag(q))                         
            u<-M_1%*%(t(Wj)%*%t(Jc) + Sigj*(Alphaj%*%Covarsj))                   ## Expected value of scores.
            Sum_Euu<-(nrow(Jc)*Sigj*M_1) + (u%*%t(u))
            
            # M-step
            Alphaj<-(u%*%t(Covarsj))%*%solve(Covarsj%*%t(Covarsj))                         ## Estimation of regression coefficients
            
            Wj<-(t(Jc)%*%t(u))%*%solve(Sum_Euu)                             ## Estimation of loadings.                
                     
            YWEu<-sum(diag(Jc%*%Wj%*%u))
            MLESigj<-(nrow(Jc)*sum(diag(Sj)) + sum(diag((t(Wj)%*%Wj)%*%Sum_Euu)) - 2*YWEu)/(p*nrow(Jc))    ## Estimation of the variance.
      		Sigj<- c(((N*p)*MLESigj + C2p)/((N*p) + Vp + 2))


            ### Observed Log Likelihood jackknife PPCCA ########
            Den<-rep(NA, nrow(J))
            Sigma<-Wj%*%t(Wj)+Sigj*diag(p)
            mumat<-Wj%*%(Alphaj%*%Covarsj) + matrix(colMeans(J), nrow=p, ncol=nrow(J), byrow=FALSE)
            for(i in 1:nrow(J))
            {
              Den[i]<-(dmvnorm(J[i,], mumat[,i], Sigma, log=TRUE))
            }
            ll[q,v]<-sum(Den)

     		#### Covergence Assessment
            converge<-Aitken(ll, lla, v, q, epsilon)
            tol<-converge[[1]]
            lla[q,v]<-converge[[2]]
         } # close while loop 

  	   ## Store jackknife parameters 
  	   Sigstore[n]<-Sigj
  	   Wstore[,,n]<-Wj  
       Alphastore[,,n]<-Alphaj 
	}# Close n loop

	## Calculate jackknife standard errors
    se_Alpha<-sqrt((N-1)^2/N*(apply(Alphastore, c(1,2), var)))
	se_W<-sqrt((N-1)^2/N*(apply(Wstore, c(1,2), var)))

	## Compute 95% confidence intervals 
    UpperCI_Alp<-Alphaopt + (qnorm(conflevel + ((1 - conflevel)/2))*se_Alpha)
    LowerCI_Alp<-Alphaopt - (qnorm(conflevel + ((1 - conflevel)/2))*se_Alpha)
	UpperCI_W<-Wopt + (qnorm(conflevel + ((1 - conflevel)/2))*se_W)
	LowerCI_W<-Wopt - (qnorm(conflevel + ((1 - conflevel)/2))*se_W)
      
    ## Report the regression coefficient CIs
    colnames(LowerCI_Alp)<-paste("LowLim_Beta_", rownames(Covars), sep="")
	colnames(UpperCI_Alp)<-paste("UpperLim_Beta_", rownames(Covars), sep="")
	tab1<-cbind(LowerCI_Alp, UpperCI_Alp)
	for(i in 1:(L+1))
	{ 
		if(i==1)
		{
			CI_Alp<-cbind(tab1[,i], tab1[,(i+L+1)])
			colnames(CI_Alp)<-colnames(tab1)[c(i, (i+L+1))]
		}else{
			tab2<-cbind(tab1[,i], tab1[,(i+L+1)])
			colnames(tab2)<-colnames(tab1)[c(i, (i+L+1))]
			CI_Alp<-cbind(CI_Alp, tab2)
		} # if
	}#i
     
    ## Confidence intervals for the first PC.
	CI_W<-cbind(LowerCI_W[,1],UpperCI_W[,1])
	colnames(CI_W)<-c("LowCI_PC1", "UppCI_PC1")

	## Determine ppm variables significantly different from zero.
	ProdciPC1<-apply(CI_W[,1:2],1,prod)
	Signifppmz<-ProdciPC1>0                       ## Selecting ppm variables with loadings on PC 1 which are significantly different from zero.
	SignifW<-as.matrix(Wopt[Signifppmz,])
	nSL<-nrow(SignifW)
	Lower<-as.matrix(LowerCI_W[Signifppmz,])
	Upper<-as.matrix(UpperCI_W[Signifppmz,])
	cat("The number of spectral bins with loadings on PC 1 significantly different from 0 is:", nSL, "\n \n")

	## Select significant ppm variables with `high' loadings 
	grid<-0.1
	cutoff<-seq(min(abs(SignifW)), max(abs(SignifW)), by=grid)

	nHL<-rep(NA,length(cutoff))
	for(l in 1:length(nHL))
	{
   		SignifHighppm<-abs(SignifW[,1])>cutoff[l]
   		SignifHighW<-SignifW[SignifHighppm,]
   		nHL[l]<-sum(SignifHighppm)
	}

	## Plot cutoff versus number of significant high loadings for user decision on value of cutoff. 
	plot(cutoff, nHL , xlab="Cutoff", ylab="Number of spectral bins selected", type="l",  cex=1.2, pch=21, col="blue", col.lab="blue", lwd=2, main="Spectral bin selection", col.main="blue")
	dat<-data.frame(Cutoff=round(cutoff,2), No.Selected=nHL)
	print(data.frame(dat))
	cat("\n \n")
	
	number<-as.numeric(ask("The user should select the number of spectral bins which have significantly high loadings. \n \n Please type the number of spectral bins you require:"))
	cutoff<-cutoff[(nHL==number)][1]
	
	SignifHighppm<-abs(SignifW[,1])>cutoff
    SignifHighW<-matrix(SignifW[SignifHighppm,], ncol=q)
	LowerH<-matrix(Lower[SignifHighppm,], ncol=q)
	UpperH<-matrix(Upper[SignifHighppm,], ncol=q)


	## Bar plot of significant high loadings
	barplot2(SignifHighW[,1], ylim=c(min(LowerH[,1],0)-1.5, max(UpperH[,1],0)+1.5), las=2, width=0.5, space=0.5, plot.grid=TRUE, ylab="PC 1 loadings", xlab="Spectral regions", names.arg = rownames(SignifW)[SignifHighppm], plot.ci = TRUE, ci.l = LowerH[,1], ci.u = UpperH[,1], font.main = 2, col="red")
	
list(q=q, sig=Sigopt, scores=Uopt, loadings=Wopt, SignifW=SignifW, SignifHighW= SignifHighW, LowerCI_W=Lower, UpperCI_W=Upper, coefficients=Alphaopt, coeffCI=CI_Alp, Cutoffs=dat, number=number, cutoff=cutoff, BIC=BIC, AIC=AIC)
} #End ppcca.metabol.jack

