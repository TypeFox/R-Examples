############################################################################################
##
##  R function for variance estimation of FMT method (Fully Moderated t-statistic)    
##  Author: Lianbo Yu, The Ohio State University                                      
##  Contact: Lianbo.Yu@osumc.edu
##  Last update: March 2012                                                             
##
############################################################################################
##
##  This function computes posterior residual variances to be used in the denominator of 
##  a moderated t-statistic from a linear model analysis of microarray data.  It is an  
##  extension of the moderated t-statistic original proposed by Smyth (2004). LOESS local 
##  regression and empirical Bayesian method are used to estimate gene specific prior
##  degrees of freedom and prior variance based on average gene intensity level. The 
##  posterior residual variance in the denominator is a weighted average of prior and 
##  residual variance and the weights are prior degrees of freedom and residual variance 
##  degrees of freedom. The degrees of freedom of the moderated t-statistic is simply
##  the sum of prior and residual variance degrees of freedom.
##
##
##  Inputs:
##
##  Amean     - vector of average log intensity levels of all genes
##  sigmasq   - vector of residual variances of all genes
##  df        - degrees of freedom for sigmasq
##  span1     - span parameter in LOESS smoothing function
##  span2     - span parameter in LOESS smoothing function
##  iter1     - iteration number in LOESS smoothing function
##  iter2     - iteration number in LOESS smoothing function
##  b         - number of genes on either side of moving average window 
##              when calculating variance of log residual variances
##
##
##  Outputs:
##   
##  df.prior  - estimated prior degrees of freedom
##  df.post   -	estimated posterior degrees of freedom
##  s2.prior  -	estimated prior variance
##  s2.post   -	estimated posterior variance
##  Ameansort - intermediate result for plotting
##  eg        - intermediate result for plotting
##  egpred    - intermediate result for plotting
##  MAvar     - intermediate result for plotting
##  tri.d0    - intermediate result for plotting
##		
##	
##
##  Example Function Call:
##
##  result <- fmt(Amean,sigmasq,df)
##
##
############################################################################################

fmt <- function(Amean, sigmasq, df, ...) UseMethod("fmt")

fmt.default <- function(Amean, sigmasq, df, span1=0.5, span2=0.95, iter1=4, iter2=4, b=20, ...) {

        library(limma)

	## LOESS smoothing of log variances
	eg <- log(sigmasq) - digamma(df/2) + log(df/2)
	egpred <- loessFit(eg, Amean, iterations=iter1, span=span1)$fitted

	## moving average calculation of variances of log variances
        N <- length(Amean)
	mat <- cbind(Amean,(eg - egpred)^2)
	order <- sort(Amean,index.return=TRUE)$ix
	matsort <- mat[order,]
	MAvar <- NULL
	for (i in 1:b) {
 	      MAvar[i] <- mean(matsort[1:(i+b),2])
	}
	for (i in (b+1):(N-b)) {
	      MAvar[i] <- mean(matsort[(i-b):(i+b),2])
	}
	for (i in (N-b+1):N) {
	      MAvar[i] <- mean(matsort[(i-b):N,2])
	}

        ## LOESS smoothing of variances of log variances
	tri.d0 <- loessFit(MAvar, matsort[,1], iterations=iter2, span=span2)$fitted - trigamma(df/2) 
	tri.d0[tri.d0<=0] <- 0.001

	## prior and posterior estimates
	df.prior.sort <- 2*trigammaInverse(tri.d0) 
        df.prior <- df.prior.sort[sort(order,index.return=TRUE)$ix]
	df.post <- df.prior + df        
   	s2.prior <- exp(egpred + digamma((df.prior)/2) - log(df.prior/2))
	s2.post <- (df.prior*s2.prior + df*sigmasq) / df.post
        
        ## output
	output <- data.frame(df=df,df.prior=df.prior,df.post=df.post,s2.prior=s2.prior,s2.post=s2.post,
                             Amean=Amean,Ameansort=matsort[,1],eg=eg,egpred=egpred,MAvar=MAvar,tri.d0=tri.d0)

        class(output) <- "fmt"
	return(output) 
}






#############################################################
##
##  R function for plotting results from FMT function      
##  Author: Lianbo Yu, The Ohio State University           
##  Contact: Lianbo.Yu@osumc.edu
##  Last update: March 2012                                 
##
##  Example Function Call:
##
##  plot(result,type="all")
##
#############################################################


plot.fmt <- function(x, type, ...) {

          if (type=="all") {
		par(mfrow = c(2,2))
		plot(x$Amean,x$df.prior,cex=0.5,col="black",xlab="Average Log Intensity",ylab="Prior Degrees of Freedom",type="p")
 		plot(x$Ameansort,x$MAvar,cex=0.5,ylab="Variance of Log Variance",xlab="Average Log Intensity",pch=1)
		points(x$Ameansort,x$tri.d0+trigamma(x$df/2),cex=0.5,col="red")
 		plot(x$Amean,x$s2.prior,cex=0.5,col="black",xlab="Average Log Intensity",ylab="Prior Variance",type="p")
		plot(x$Amean,x$eg,cex=0.5,xlab="Average Log Intensity",ylab="Log Variance",pch=1)
 		points(x$Amean,x$egpred,col="red",cex=0.5)
          }
          if (type=="priordf") {
		par(mfrow = c(1,1))
		plot(x$Amean,x$df.prior,cex=0.5,col="black",xlab="Average Log Intensity",ylab="Prior Degrees of Freedom",type="p")
          }
          if (type=="varoflogvar") {
		par(mfrow = c(1,1))
 		plot(x$Ameansort,x$MAvar,cex=0.5,ylab="Variance of Log Variance",xlab="Average Log Intensity",pch=1)
		points(x$Ameansort,x$tri.d0+trigamma(x$df/2),cex=0.5,col="red")
          }
          if (type=="priorvar") {
		par(mfrow = c(1,1))
 		plot(x$Amean,x$s2.prior,cex=0.5,col="black",xlab="Average Log Intensity",ylab="Prior Variance",type="p")
          }
          if (type=="logvar") {
		par(mfrow = c(1,1))
		plot(x$Amean,x$eg,cex=0.5,xlab="Average Log Intensity",ylab="Log Variance",pch=1)
 		points(x$Amean,x$egpred,col="red",cex=0.5)
          }
}



