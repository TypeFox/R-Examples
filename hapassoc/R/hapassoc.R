# Filename: hapassoc.R
# Version : $Id: hapassoc.R,v 1.60 2012/04/05 04:09:37 mcneney Exp $

# hapassoc- Inference of trait-haplotype associations in the presence of uncertain phase
# Copyright (C) 2003  K.Burkett, B.McNeney, J.Graham

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

########################################################################

hapassoc<-function(form, haplos.list, baseline="missing", family=binomial(), 
                   design="cohort", disease.prob=NULL, freq=NULL, maxit=50, 
                   tol=0.001, start=NULL, verbose=FALSE){
environment(form)<-environment()#set envir of formula to envir w/i hapassoc function
call <-  match.call()

if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
if (is.function(family))
    family <- family()
if(design=="cc" && family$family!="binomial") {
  warning(paste("design is case-control, but family is ",family$family,"; using binomial instead",sep=""))
  family<-binomial()
}
#if(family$family=="binomial") family<-binomialNoWarn() 
   # so we don't issue warnings in M step

haplos.names<-names(haplos.list$initFreq)

 # Initial freq values, if no freq specified use initFreq
 if (!is.null(freq)) {
   names(freq)<-haplos.names
 } else {
   freq<-haplos.list$initFreq
 }

 if(baseline=="missing") {
   #no baseline specified
   baseline<-haplos.names[which.max(freq)] 
 }

 if ( length( setdiff( colnames(haplos.list$haploDM), all.vars(form)[2:length(all.vars(form))]))==0){
   #Too many haplotypes specified in model formula. Must determine a baseline category and 
   #create a new model formula. 
   column.subsethaplo <- colnames(haplos.list$haploDM) != baseline
   resp <- as.character(all.vars(form)[1])
   column.subsetnonhaplo <- intersect(colnames(haplos.list$nonHaploDM),all.vars(form)[2:length(all.vars(form))])
   form.rhs <- paste(c(names(haplos.list$nonHaploDM[, column.subsetnonhaplo, 
            drop = FALSE]), names(haplos.list$haploDM)[column.subsethaplo]), 
            collapse = "+")
   form <- formula(paste(resp, "~", form.rhs))
   cat(paste("Formula contained all possible haplotypes; using ",baseline," as reference 
category\n\n"))
 }

 if(!is.na(all.vars(form)[2]) && all.vars(form)[2]==".") {
   #Then formula is of form "y ~ ." 
   column.subsethaplo <- colnames(haplos.list$haploDM)!=baseline
   resp<-as.character(all.vars(form)[1])
   column.subsetnonhaplo <- colnames(haplos.list$nonHaploDM)!=resp
   form.rhs<-paste(c(names(haplos.list$nonHaploDM[,column.subsetnonhaplo,drop=FALSE]),
				   names(haplos.list$haploDM)[column.subsethaplo]),collapse="+")
   form<-formula(paste(resp,"~",form.rhs))
 }
 hdat <- cbind(haplos.list$nonHaploDM, haplos.list$haploDM)
 colnames(hdat)<- c(colnames(haplos.list$nonHaploDM),
                    colnames(haplos.list$haploDM))
 ID <- haplos.list$ID
 N<-round(sum(haplos.list$wt))
 wts<-haplos.list$wt
 
 # Get the haplotype columns
 haplos<-haplos.list$haploDM
 haploMat<-haplos.list$haploMat
 allHaps<-c(haploMat[,1],haploMat[,2]) #Needed later in hapassoc loop for wt calcs
 # Initial beta values calculated from augmented data set
 if (family$family=="binomial"){
	 regr=tryCatch( glm(form, family=family, 
					 data=hdat, weights=wts, start=start),
			 
			 warning = function(w) {
				 
				 teststring="non-integer #successes in a binomial glm!"
				 if (teststring==conditionMessage(w)){
					 suppressWarnings(glm(form, family=family, 
									 data=hdat, weights=wts, start=start))
				 } else {
					 glm(form, family=family, data=hdat, weights=wts, start=start)
				 }
			 }
	 
	 )
 } else { 
	 regr<-glm(form, family=family, data=hdat, weights=wts, start=start)
 }
 
 response<-regr$y #Change added Nov.2/04 to extract response from fitted model
 beta<-regr$coef
 fits<-regr$fitted.values
 betadiff<-1
 it<-1
 num.prob<-vector(length=nrow(hdat))
 
 # The hapassoc loop for case-control data

 if (design=="cc")    ## use the hybrid method for case-control study
 {  		
     n1 <- sum(response[!duplicated(ID)])
     n0 <- N - n1

     nhaps <- length(freq) 	
     nhaps.DM <- length(haplos.list$haploDM) 	

   # Construct the design matrix "Des.Mat" for the psuedo-sample.     
     mf <- model.frame(form, hdat)
     Des.Mat <- model.matrix(form, mf, contrasts)  

   # Construct the orginal non-haplotype data "nonHap" for the subjects.     
   # The non-haplo data may include factors, but we require the numeric
   # columns of the design matrix. We will extract these from Des.Mat by
   # trimming off haplotype columns. Haplotype columns appear last in Des.Mat
     nonHap <- haplos.list$nonHaploDM[!duplicated(ID),,drop=FALSE]
     #nonHap<- as.matrix(nonHap)     

   # Construct the design matrices (array) "DesMat.Arr" for numerator of the 
   # effective sample-sizes (see formula (13) of Spinka et al.).
     nonHap.Mat <- matrix(nrow=nrow(nonHap)*length(freq), ncol=ncol(nonHap))
     nonHap.Mat<-data.frame(nonHap.Mat)
     for (j in 1:ncol(nonHap)) {
        if(is.numeric(nonHap[,j])) {
          nonHap.Mat[,j] <- c(outer(rep.int(1,length(freq)), nonHap[,j]))
        } else { #assume factor
          if(!is.factor(nonHap[,j])) stop("non-numeric, non-factor covariate")
          n.nonHapj<-as.numeric(nonHap[,j])
          nonHap.Matj <- c(outer(rep.int(1,length(freq)), n.nonHapj))
          nonHap.Mat[,j]<-factor(nonHap.Matj,labels=levels(nonHap[,j]))
        }
     }
     DesMat.Arr <- array(data=NA, dim=c(nrow(nonHap.Mat), length(beta), nhaps))
     for (j in 1:nhaps)
     {
#       Hap.Mat <- matrix(rep.int(diag(1, nhaps), N), ncol=nhaps, byrow=TRUE)
       Hap.Mat <- matrix(rep.int(diag(1, length(freq)), N), ncol=length(freq), byrow=TRUE)
       Hap.Mat[,j] <- Hap.Mat[,j] + 1
       indPooled <-  match(haplos.list$pooled,names(freq))  
       if(!any(is.na(indPooled))) { #then there is pooling to do
         pooled <- vector("numeric",nrow(Hap.Mat))
         for (k in 1:length(indPooled))
           pooled <- pooled + Hap.Mat[,indPooled[k]]
         dm <- data.frame(nonHap.Mat, Hap.Mat[,-indPooled], pooled)
       }
       else  dm <- data.frame(nonHap.Mat, Hap.Mat)
       dimnames(dm)[[2]] <-c(dimnames(haplos.list$nonHaploDM)[[2]], 
                             colnames(haplos.list$haploDM) )
       mf <- model.frame(form, dm)
       DesMat.Arr[,,j] <- model.matrix(form, mf, contrasts)
     }
     DMA.ID <- c(outer(rep.int(1,nhaps),unique(ID)))

   # Construct the design matrix "Denom.Mat" for the denominator of the 
   # effective sample-sizes  

     DiploMat <- matrix(nrow=length(freq)*(length(freq)+1)/2, ncol=length(freq))
     myInd <- 1
     for (j in 1:length(freq))  {
       diag.mat <- diag(1,length(freq))
       diag.mat[,j] <- diag.mat[,j] + 1
       DiploMat[myInd:(myInd + length(freq)-j),] <- diag.mat[j:length(freq),]
       myInd <- myInd + length(freq)-j+1
     }
     if(!any(is.na(indPooled))) {
       pooled <- vector("numeric",nrow(DiploMat))
       for (k in 1:length(indPooled))
         pooled <- pooled + DiploMat[,indPooled[k]]
       DiploMat <- data.frame(DiploMat[,-indPooled], pooled)
     }

     Denom.Mat<-matrix(NA,nrow=nrow(DiploMat)*nrow(nonHap),ncol=ncol(nonHap))
     Denom.Mat<-as.data.frame(Denom.Mat)
     for (j in 1:ncol(nonHap)) {
       if(is.numeric(nonHap[,j])) {
         Denom.Mat[,j] <- c(outer(rep.int(1,nrow(DiploMat)), nonHap[,j]))  
       } else {
         n.nonHapj<-as.numeric(nonHap[,j])
         Denom.Matj<- c(outer(rep.int(1,nrow(DiploMat)), n.nonHapj))  
         Denom.Mat[,j]<-factor(Denom.Matj,labels=levels(nonHap[,j]))
       }
     }
     for (j in 1:nhaps.DM)
       Denom.Mat <- data.frame(Denom.Mat, rep.int(DiploMat[,j], N))
       
     dimnames(Denom.Mat)[[2]] <- dimnames(dm)[[2]]    
     mf <- model.frame(form, Denom.Mat)
     Denom.Mat <- model.matrix(form, mf, contrasts)
     DenomMat.ID <- c(outer(rep.int(1, nrow(DiploMat)), unique(ID)))
     
 		
     while ( (it<maxit) && (betadiff>tol) ){
 	 	
    # Multiplicative const for haplo probs: 1 for homo, 2 for het

      haplo.probs<-rep.int(1,nrow(haplos))+isMultiHetero(haplos.list)
      haplo.probs <- haplo.probs*freq[haploMat[,1]]*freq[haploMat[,2]]
      phi<-mlPhi(regr) #Compute ML estimate of dispersion param
      if(is.null(phi)) { #no converergence in ml estimate of phi
           break() #hapassoc will throw a warning of non-convergence
      }

      num.prob <- exp(response*(Des.Mat%*%beta))
      if (!is.null(disease.prob))  # if Pr(D=1) is known, beta0 can be calulated 
      {   
         beta0 <- beta[1] - log(n1/n0) + log(disease.prob/(1-disease.prob))
         num.prob <- num.prob/(1+exp(1+Des.Mat%*%c(beta0,beta[-1]))) 
      }
      num.prob <- as.vector(num.prob)*haplo.probs   

    # E step: Calculate the weights for everyone
    # Use the ID to determine the number of pseudo-individuals in the 
    # denominator probability
      wts<-vector("double", length(wts))
      wts<-.C("getWts", as.integer(nrow(hdat)), as.integer(ID),
              wts = as.double(wts), as.double(as.vector(num.prob)),
              PACKAGE="hapassoc")
      wts <- wts$wts

     
     if (min(freq)<1.0e-8 | min(wts)<1.0e-8)  
     {  # if some haplotype frequency is estimated to be 0, 
     	 # throw a non-convergence warning
       break()
     }
    
    # M step: Find new estimates using GLM and weighted haplotype counts

    # Find the new betas using old betas as starting value
    
	if (family$family=="binomial"){
		regr=tryCatch( glm(form, family=family, data=hdat, weights=wts, 
						control=glm.control(epsilon=1e-08),start=beta),
					 
					 warning = function(w) {
						 
						 teststring="non-integer #successes in a binomial glm!"
						 if (teststring==conditionMessage(w)){
							 suppressWarnings(glm(form, family=family, data=hdat, weights=wts, 
											 control=glm.control(epsilon=1e-08),start=beta))
						 } else {
							 glm(form, family=family, data=hdat, weights=wts, 
									 control=glm.control(epsilon=1e-08),start=beta)
						 }
					 }
			 
			 )
		 } else { 
			 regr <- glm(form, family=family, data=hdat, weights=wts,
					 control=glm.control(epsilon=1e-08),start=beta)
		 }		 
		 
		 
		 
      betaNew<-regr$coef
      fits<-regr$fitted.values
      betadiff<-max( max(abs(beta-betaNew), na.rm=TRUE),
	       max(abs(beta-betaNew)/(0.1+abs(betaNew)), na.rm=TRUE) )
      beta<-betaNew

    # Calculate the expected haplotype counts "Numer.Freq" 
    # in formula (13) of Spinka el al. (2005) 
      Numer.Freq <- tapply(c(wts,wts),allHaps,sum)

    # Calculate the effective sample-sizes "Denom.Freq" 
    # in formula (13) of Spinka el al. (2005)
      DenomFreq.Nu <- NULL
      for (j in 1:nhaps)
      {
        rr <- r.Omega(DesMat.Arr[,,j], beta, ncases=n1, ncontrols=n0, disease.prob) 
        rr <- 2*rep.int(freq,N)*rr 

        factors <- unique(DMA.ID)
        sums <- vector("double",length(factors))
        rr <- .C("tapply_sum", as.integer(length(rr)), as.integer(factors), 
           as.double(rr), as.integer(DMA.ID), sums = as.double(sums),
           PACKAGE="hapassoc")
        rr <- rr$sums
        names(rr) <- factors

        DenomFreq.Nu <- cbind(DenomFreq.Nu, rr)
      }
      DenomFreq.Nu <- data.frame(DenomFreq.Nu)
      dimnames(DenomFreq.Nu)[[2]] <- haplos.names
      DenomFreq.De <- rep.int(get.diplofreq(freq),N)*r.Omega(Denom.Mat,beta, 
                      ncases=n1,ncontrols=n0,disease.prob)


      factors <- unique(DenomMat.ID)
      sums <- vector("double",length(factors))
      DenomFreq.De <- .C("tapply_sum", as.integer(length(DenomFreq.De)), 
                       as.integer(factors), as.double(DenomFreq.De), 
                       as.integer(DenomMat.ID), sums = as.double(sums),
                       PACKAGE="hapassoc")
      DenomFreq.De <- DenomFreq.De$sums      

      Denom.Freq <- apply(DenomFreq.Nu/DenomFreq.De, 2, sum)  
    
    # Update the haplotype frequencies  
      freq <- Numer.Freq/Denom.Freq
      freq <- freq/sum(freq)
    
      if (verbose) 
	      cat("iteration",it,": value of convergence criterion =",betadiff,"\n")

      it<-it+1
   } # end of while
 } # end of if (design=="cc")

 else   # If design is not "cc", use original hapassoc code.
 { 
   while ( (it<maxit) && (betadiff>tol) ){
   
        # Vector of P(Y|X)*P(X) probability 
	# If the person is unaffected, P(Y=0) is 1-fitted value
	# otherwise it is the fitted value 
   
        # Multiplicative const for haplo probs: 1 for homo, 2 for het
        haplo.probs<-rep.int(1,nrow(haplos))+isMultiHetero(haplos.list)
        haplo.probs <- haplo.probs*freq[haploMat[,1]]*freq[haploMat[,2]]

	phi<-mlPhi(regr) #Compute ML estimate of dispersion param
        if(is.null(phi)) { #no converergence in ml estimate of phi
           break() #hapassoc will throw a warning of non-convergence
        }
        num.prob<-pYgivenX(response,fits,phi,family)*haplo.probs

	# E step: Calculate the weights for everyone
	# Use the ID to determine the number of pseudo-individuals in the 
	# denominator probability

        wts<-vector("double", length(wts))
        wts<-.C("getWts", as.integer(nrow(hdat)), as.integer(ID),
                 wts = as.double(wts), as.double(as.vector(num.prob)),
                 PACKAGE="hapassoc")
        wts <- wts$wts

	# M step: Find new estimates using GLM and weighted haplotype counts

	# Find the new betas using old betas as starting value
        
	if (family$family=="binomial"){
		regr=tryCatch( glm(form, family=family, data=hdat, weights=wts, 
						control=glm.control(epsilon=1e-08),start=beta),
						
						warning = function(w) {
							
							teststring="non-integer #successes in a binomial glm!"
							if (teststring==conditionMessage(w)){
								suppressWarnings(glm(form, family=family, data=hdat, weights=wts, 
												control=glm.control(epsilon=1e-08),start=beta))
							} else {
								glm(form, family=family, data=hdat, weights=wts, 
										control=glm.control(epsilon=1e-08),start=beta)
							}
						}
				
				)
	} else { 
		regr <- glm(form, family=family, data=hdat, weights=wts,
				control=glm.control(epsilon=1e-08),start=beta)
	}	
			
    betaNew<-regr$coef
   	fits<-regr$fitted.values
	betadiff<-max( max(abs(beta-betaNew), na.rm=TRUE),
	max(abs(beta-betaNew)/(0.1+abs(betaNew)), na.rm=TRUE) )   	
	beta<-betaNew

	# Find the new freqs, weighted sum of haplotypes
    freq <- tapply(c(wts,wts),allHaps,sum)/(2*N)

	if(verbose) 
	  cat("iteration",it,": value of convergence criterion =",betadiff,"\n")

    it<-it+1
   } # end of while
 } # end of else: design != "cc" 
 
 	 
 if(betadiff>tol) { #did not converge
    warning(paste("No convergence in hapassoc in",it,"iterations\n")) 
    ans<-list(converged=FALSE)
    class(ans)<-"hapassoc"
    return(ans)
 }

 #freq is currently an array which causes problems in the calculations below
 freq <- as.matrix(freq) 

 # Added more elements to EMresults list for later call to
 # log.likelihood and renamed dispersion (KB 18/02/2006) and
 # freq (KB 28/07/2008) 
 EMresults <- list(beta=beta, freq=freq, fits=fits, wts=wts, ID=ID,
                   glm.final.fit=regr,dispersion=phi, family=family,
                   response=response)


 # Compute the log-likelihood so that results can be returned
 # from function (KB 18/02/2006)
 if(design!="cc") {
   loglik <- loglikelihood(haplos.list,EMresults)
 } else {
   loglik<-NA
 }

 names(haplos.list)[names(haplos.list) == "freq"] <- "gamma"
 names(EMresults)[names(EMresults) == "freq"] <- "gamma"
 var.est <- EMvar(haplos.list, EMresults, family)

 # Added the model equation and log-likelihood as elements
 # returned from function (KB 18/02/2006)
 ans<-list(it=it, beta=beta, freq=freq, fits=fits, wts=wts, ID=ID,
           var=var.est, dispersion=phi, family=family, response=response,
           converged=TRUE, model=form,loglik=loglik, call=call)

 class(ans)<-"hapassoc"
 return(ans)

}

## Other functions called in hapassoc:

#############################################################################
## Function to calculate r.Omega
r.Omega <- function(des.mat, beta, ncases, ncontrols, disease.prob=NULL)
{
 # Assuming a rare disease 	
   rr <- 1 + exp(des.mat%*%beta)  
 # If Pr(D=1) is known, beta_0 can be calculated
   if (!is.null(disease.prob))    
   {
 	  beta0 <- beta[1] - log(ncases/ncontrols) + log(disease.prob/(1-disease.prob))
 	  rr <- rr / (1 + exp(des.mat%*%c(beta0,beta[-1])))
   }	
   return(as.vector(rr))
}
 
#############################################################################
# Function to calcuate the diplotype frequencies, given haplotype frequencies 
get.diplofreq <- function(freq)
{
  nhaplos <- length(freq)
  dipprob <- 2*kronecker(freq,t(freq))
  diag(dipprob) <- diag(dipprob)/2
  dipprob.vec<-NULL
  for(i in 1:nhaplos) 
    dipprob.vec<-c(dipprob.vec,dipprob[i:nhaplos,i])
  return (dipprob.vec)
}


########################################################################

mlPhi<-function(myfit,maxit=30,tol=1e-5) {
  if(myfit$family$family=="binomial" || myfit$family$family=="poisson") {
    return(1) #dispersion set to one
  }
  if(myfit$family$family=="gaussian") {
    return(summary(myfit)$deviance/sum(myfit$prior.weights))
  }
  if(myfit$family$family=="Gamma") { #Need to do Newton-Raphson
    return(mlPhiGamma(myfit,maxit,tol))
  }
  #Else, we don't support the family passed in
  stop(paste("ML estimation of dispersion not supported for family",
                myfit$family$family))
}

########################################################################

mlPhiGamma<-function(myfit,maxit,tol) {
  
  # Need to do Newton-Raphson, use moment estimate of phi to get
  # starting value for N-R 

  dev<-myfit$deviance
  n<-sum(myfit$prior.weights)
  phiMoment<-summary(myfit)$dispersion #moment estimate of phi
  ## Looking for root of score equation as a function of x=1/phi
  ## Function f is the score equation and fp is its derivative
  f<-function(x,n,dev) { return(dev+2*n*(digamma(x)-log(x))) }
  fp<-function(x,n) { return(2*n*(trigamma(x)-1/x)) }
  diff<-1; i<-1  # initialization
  xold<-1/phiMoment #starting value for N-R 
  while(i<=maxit && diff>tol) {
    xnew<-xold - f(xold,n,dev)/fp(xold,n)
    if(is.na(xnew)) {
       warning("No convergence for ML estimate of Gamma scale parameter")
       return(NULL)
    }
    diff<-abs(xnew-xold); xold<-xnew; i<-i+1
  }
  if(i>maxit) 
    warning("No convergence for ML estimate of Gamma scale parameter")
  return(1/xnew)
}

########################################################################

pYgivenX<-function(y,mu,phi,family){

  #Calculate P(Y|X) to be used in the weights. We use P(Y|X) in 
  #expressions like P(Y_i|X_i^j)P(X_I^j)/sum_k{ P(Y_i|X_i^k)P(X_i^k) }
  #so factors that are the same for all X_i^j can be ignored.
  #The deviance resids can be used to get P(Y|X) up to constants:
  #Assuming weights of 1 
  #dev.resid = 2*phi*(l(y,phi) - l(mu,phi)) where l(mu,phi) is the log-lik
  #using mean mu and l(y,phi) is the log-likelihood in the saturated
  #model using y as the mean.
  #So   exp(-dev.resid/(2*phi)) = P(y|x)*exp(-l(y,phi))
  #Call the dev.resid function with a weight of 1 

  return(exp(-1*family$dev.resid(y,mu,wt=1)/(2*phi)))
} 


## EMvar Functions
 
########################################################################

EMvar<-function(haplos.list, results, family)  {
 # Get the results of the last iteration
 beta<-results$beta[!is.na(results$beta)]
 betanames<-names(beta)
 num.beta<-length(beta)

 gammanames<-rownames(results$gamma)
 gamma<-as.vector(results$gamma)
 num.gamma<-length(gamma)

 weights<-results$wt
 ID <- haplos.list$ID

 final.regr <- results$glm.final.fit

 # Set up the Y vectors and P vector of fitted values
 # one for the augmented data set, one for just the missing
 
 missing<-(weights<1)

 yi.full<-results$response
 fits.full<-results$fits
 
 yi.mis<-yi.full[missing]
 fits.mis<-fits.full[missing]

 # Set up the matrices needed to find the variance:

   ## Design matrices

   Xreg <- model.matrix(final.regr)
   colnames(Xreg)<-betanames
   Xreg.mis<-Xreg[missing,]

   ## Haplotype matrix
   #Instead of assuming the glm model is additive in the haplotypes, 
   #make the Xgen matrix from scratch here. Used to be
   #   Xgen<-as.matrix(haplos.list$haploDM)
   haploMat<-haplos.list$haploMat
   #each application of outer in next line returns a matrix of T/F with
   # (i,j)th element true if haploMat[i,k]==gammanames[j]; k=1,2
   #Adding these T/F matrices coerces to 1/0 s first before adding
   Xgen<-outer(haploMat[,1],gammanames,"==")+outer(haploMat[,2],gammanames,"==")

   Xgen.mis <- as.matrix(Xgen[missing,])
 
   ## Weight matrices

   W.full<-diag(weights)
   W.mis<-diag(weights[missing])

   ## likelihood derivative matrices

   H.mis<-diag(yi.mis-fits.mis)/results$dispersion

   ## Haplotype frequency and derivative matrices

   ones.mat<-matrix(1,nrow=(num.gamma-1), ncol=(num.gamma-1))
   num.haplo<-vector(length=num.gamma)
   for (i in 1:num.gamma){
	num.haplo[i]<-sum(weights*Xgen[,i])}
   der.gamma<-1/gamma
   der2.gamma<-1/gamma^2

   G1<-diag(der.gamma)
   G<-G1[,1:(num.gamma-1),drop=FALSE]
   G[num.gamma,]<-G1[num.gamma, num.gamma]

   ## Block diagonal matrix

   ID.mis <- ID[missing]
   i <- 1
   B<-matrix(0,nrow=sum(missing), ncol=sum(missing) )

   while (i<length(ID.mis)){
	pseudo.index<-ID.mis==ID.mis[i]
	ones.vec<-rep.int(1, sum(pseudo.index)*sum(pseudo.index))
	block<-matrix(ones.vec, nrow=sum(pseudo.index), ncol=sum(pseudo.index))
	block.size<-sum(pseudo.index)
	B[i:(i+block.size-1), i:(i+block.size-1)]<-block
	i<-i+sum(pseudo.index)
   }

 # Calculate the complete data expected information (Ic)

   ## Top block

   Ic.reg<-solve(summary(final.regr)$cov.unscaled)/results$dispersion

   ## Middle block if dispersion was estimated -- Ic.phi requires Sphi

   Sphi<-SPhi(final.regr,results$dispersion) #will be NULL if phi not est'd
   if(!is.null(Sphi)) { 
     Ic.phi<-IPhi(final.regr,Sphi,results$dispersion)
   } else { 
     Ic.phi<-NULL 
   }

   ## Bottom block

   if(num.gamma>2) {
     Ic.cov<- diag(num.haplo[1:(num.gamma-1)]*der2.gamma[1:(num.gamma-1)])+der2.gamma[num.gamma]*num.haplo[num.gamma]*ones.mat
   } else { # first term in above sum is a scalar s, 
            # and diag(s) returns an s-by-s matrix -- not what we want
     Ic.cov<- num.haplo[1:(num.gamma-1)]*der2.gamma[1:(num.gamma-1)]+der2.gamma[num.gamma]*num.haplo[num.gamma]*ones.mat
   }

   ## Full block diagonal matrix (combine top, middle and bottom blocks)

   if(!is.null(Ic.phi)) {
     Ic<-matrix(0,nrow=(nrow(Ic.reg)+1+nrow(Ic.cov)),
                ncol=(nrow(Ic.reg)+1+nrow(Ic.cov))) #+1 is for phi
     Ic[1:num.beta,1:num.beta]<-Ic.reg
     Ic[(num.beta+1),(num.beta+1)]<-Ic.phi
     Ic[(num.beta+2):(num.beta+num.gamma),
        (num.beta+2):(num.beta+num.gamma)]<-Ic.cov
   } else {
     Ic<-matrix(0,nrow=(nrow(Ic.reg)+nrow(Ic.cov)),
                ncol=(nrow(Ic.reg)+nrow(Ic.cov)))
     Ic[1:num.beta,1:num.beta]<-Ic.reg
     Ic[(num.beta+1):(num.beta+num.gamma-1),
        (num.beta+1):(num.beta+num.gamma-1)]<-Ic.cov
   }
     

 # Calculate the missing data information (Imis), only over missing

   ## S_beta
 
   Sbeta<-H.mis%*%Xreg.mis

   ## S_phi

   if(!is.null(Sphi)) {
     Sphi<-Sphi[missing]
   }
 
   ## S_gamma
 
   Sgamma<-Xgen.mis%*%G

   ## Score matrix

   if(!is.null(Sphi)) # a fix to bug submitted by Kelly on Feb 8, 2006 -s.b.
   S<-cbind(Sbeta,Sphi,Sgamma)
   else
   S<-cbind(Sbeta,Sgamma)

   ## Calculate Imis from the other matrices

   Imis<-t(S)%*%(W.mis-W.mis%*%B%*%W.mis)%*%S  

 # Calculate the information and var/cov matrix

 Info<-Ic-Imis
 varEM<-solve(Info)

 gammanames <- paste("f", gammanames, sep="")
 if(!is.null(Sphi)) {
   colnames(varEM)<-c(betanames, "phi", gammanames[1:(num.gamma-1)])
 }else {
   colnames(varEM)<-c(betanames, gammanames[1:(num.gamma-1)])
 }
 rownames(varEM)<-colnames(varEM)


 return(varEM)

}


## Other functions called in EMvar:

########################################################################

SPhi<-function(myfit,phi) {
  if(myfit$family$family=="binomial" || myfit$family$family=="poisson") {
    return(NULL) #dispersion set to one so score not defined
  }
  if(myfit$family$family=="gaussian") {
    return(SPhiGaussian(myfit,phi))
  }
  if(myfit$family$family=="Gamma") { 
    return(SPhiGamma(myfit,phi))
  }
}

########################################################################

SPhiGaussian<-function(myfit,phi) {
  return( (myfit$y-myfit$fitted)^2/(2*phi^2) - 1/(2*phi) )
}

########################################################################

SPhiGamma<-function(myfit,phi) {
  x<-1/phi
  y<-myfit$y
  mu<-myfit$fitted.values
  return( (digamma(x)-log(x)+(y-mu)/mu-log(y/mu))*x^2 )
}

########################################################################

IPhi<-function(myfit,score,phi) {
  if(myfit$family$family=="binomial" || myfit$family$family=="poisson") {
    return(NULL) #dispersion set to one so phi not estimated
  }
  if(myfit$family$family=="gaussian") {
    return(IPhiGaussian(myfit,score,phi))
  }
  if(myfit$family$family=="Gamma") {
    return(IPhiGamma(myfit,score,phi))
  }
}

########################################################################

IPhiGaussian<-function(myfit,score,phi) {
  wts<-myfit$prior.weights
  n<-sum(wts)
  return( 2*(t(score)%*%wts)/phi + n/(2*phi^2) )
}

########################################################################

IPhiGamma<-function(myfit,score,phi) {
  x<-1/phi
  wts<-myfit$prior.weights
  n<-sum(wts)
  return( 2*x*(t(score)%*%wts) + n*x^4*(trigamma(x) - 1/x) )
}

