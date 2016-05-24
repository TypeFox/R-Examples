ebBHF <- function(formula,dom,selectdom,Xnonsample,MC=100,data,transform="BoxCox",lambda=0,constant=0,indicator)
{
   result <- list(eb=NA, fit=NA)
   
   if (missing(indicator))
      stop("Argument indicator must be included: indicator=function")

   if (transform!="BoxCox" & transform!="power")
      stop("Argument transform= \"",transform, "\" must be \"BoxCox\" or \"power\".")

   if (!missing(data))
   {
      formuladataall <- model.frame(formula,na.action = na.pass,data)
      formuladata <- model.frame(formula,na.action = na.omit,data)
      Xs  <- model.matrix(formula,data)        
      dom <- data[,deparse(substitute(dom))]
   } else
   {
      formuladataall <- model.frame(formula,na.action = na.pass)
      formuladata <- model.frame(formula,na.action = na.omit)
      Xs <- model.matrix(formula)        

   }
   if (is.factor(dom))
      dom <- as.vector(dom)

   if (nrow(formuladataall)!=length(dom))
     stop("   Arguments formula [rows=",nrow(formuladataall),"] and dom [rows=",length(dom),"] must be the same length.\n")


   welfare <- formuladata[,1]
   intercept <- (attributes(Xs)$assign[1]==0) 


   # Delete rows with NA values
   omitted     <- na.action(formuladata)
   nomitted    <- length(omitted)
   rowNAdomini <- which(is.na(dom))
   if (nomitted>0)
      dom  <- dom[-omitted]

   rowNAdom <- which(is.na(dom))
   nrowNAdom <- length(rowNAdom)
   if (nrowNAdom>0)
   {
      welfare <- welfare[-rowNAdom]
      Xs      <- Xs[-rowNAdom,]
      dom     <- dom[-rowNAdom]
   }


   # Domains for which EB estimators of indicator are required (target domains)
   if (missing(selectdom))
      selectdom <- unique(dom)
   else
   {
      if (is.factor(selectdom))
         selectdom <- as.vector(selectdom)
      selectdom <- unique(selectdom)
   }

   I    <- length(selectdom)
   n    <- length(welfare)      # Total number of sample observations
   udom <- unique(dom)
   D    <- length(udom)         # Total number of domains with sample data


   ## Test domain codes in selectdom are defined in Xnonsample
   noselectdom_Xnonsample <- NULL
   selectdomfinded        <- NULL
   udomXnonsample <- unique(Xnonsample[,1])
   for (i in 1:I)
   {
      finded <- any(selectdom[i]==udomXnonsample)
      if (finded==TRUE)
         selectdomfinded <- c(selectdomfinded,selectdom[i])
      else
         noselectdom_Xnonsample <- c(noselectdom_Xnonsample,selectdom[i])
   }
   ## warning message
   length_noselectdomXnonsample <- length(noselectdom_Xnonsample)
   if (length_noselectdomXnonsample>0)
   {
      text <- "The following domain codes (selectdom) are not defined in population \n   domain codes (Xnonsample[,1]).\n     "
      noselectdom_Xnonsample <- sort(noselectdom_Xnonsample)
      for (i in 1:length_noselectdomXnonsample)
            text <- paste(text,noselectdom_Xnonsample[i])
      warning(text)
      selectdom_ini <- selectdom
      selectdom <- selectdomfinded
      I <- length(selectdom)
   }

# A list containing a number of matrices equal to the length of selectdom.
# Matrix \code{i} must contain in each column the values of each of 
# \code{p} auxiliary variables for the out-of-sample units in \code{i}-th 
# selected domain

   p <-dim(Xs)[2]           # Number of auxiliary variables (included intercept if it is used)
   rd<-rep(0,I)
   nd<-rep(0,I)                   

   # Transform welfare after adding a constant to make it positive
   #m <- max(c(0,(-1)*min(welfare)+1))  
   m <- constant

   welfaret <- welfare+m
   checknegatives <- sum(welfaret<0)
   if (checknegatives>0)
   {
#      cat("Negativos",checknegatives,"\n")
#      cat("Minimo:",min(welfare),"\n")
      stop("   Argument constant=",constant," must be greater than ", (-1)*min(welfare))
   }
   ys       <- bxcx(welfaret,lambda=lambda,type=transform)

   # Fit the nested-error model to sample data by REML method using function lme from library nlme. 
   fit.EB<-lme(ys~-1+Xs,random=~1|as.factor(dom),method="REML")

   # Save some of the results of the fitting method in variables
   betaest   <-matrix(fixed.effects(fit.EB),nrow=p,ncol=1) # Vector of model coefficients (size p)
   upred     <-random.effects(fit.EB)                  # Predicted random effects: Watch out! It is not a vector, it is a matrix with 1 column
   sigmae2est<-fit.EB$sigma^2                          # Estimated error variance
   sigmau2est<-as.numeric(VarCorr(fit.EB)[1,1])        # VarCorr(fit2) is the estimated cov. matrix of the model random components
   sqrtsigmae2est<-sqrt(sigmae2est)

   # Create a list object containing different results from the model fit.
   Resultsfit<-list(summary=summary(fit.EB),fixed=fixed.effects(fit.EB),
                    random=upred,
                    errorvar=sigmae2est,
                    refvar=sigmau2est,
                    loglike=fit.EB$logLik,
                    residuals=fit.EB$residuals[1:n])


   # EB method starts: Generate MC vectors of non-sample values of the response 
   # from their conditional distribution given the sample data and calculate 
   # empirical values of EB estimators.

   ebestimator <- matrix(0,nrow=I,ncol=MC)  # Matrix with results for the MC simulations in the EB method
   ebestimator.EB <-rep(0,I)                # Vector with final EB estimators

   for (i in 1:I)      # Cycle for target domains
   {    
      # Matrix with the values of the p auxiliary variables for the out-of-sample 
      # observations in the i-th target domain

      d     <- selectdom[i]    # target domain code
      Xrd<-as.matrix(Xnonsample[(Xnonsample[,1]==d),-1])
      if (intercept==TRUE)
         Xrd <- cbind(1,Xrd)

      nd[i]<-sum(dom==d)
      rd[i]<-dim(Xrd)[1]

      posd  <- which(udom==d)  
      if (length(posd)!=0){
         ysd    <- ys[dom==d]          # Get sample values for target domain
         mudpred<-Xrd%*%betaest+upred[as.character(d),1]       # Compute conditional means for out-of-sample units

         # The conditional distribution of (non-sample data given sample data) in the 
         # EB method can be expressed as a new nested-error model with different random 
         # effects variance. We calculate this random effects variance (called sigmav2)
         gammad<-sigmau2est/(sigmau2est+sigmae2est/nd[i])

      } else
      {
         ysd    <- NULL               # Get sample values for target domain
         mudpred<-Xrd%*%betaest       # Compute conditional means for out-of-sample units
         gammad <- 0
      }
      sigmav2<-sigmau2est*(1-gammad)
      sqrtsigmav2<-sqrt(sigmav2)
      for (ell in 1:MC)   ### Start of Monte Carlo simulations for EB method
      {   
         vd<-rnorm(1,0,sqrtsigmav2)        # Generate random effect for target domain d
         ed<-rnorm(rd[i],0,sqrtsigmae2est) # Generate random errors for all out-of-sample units in target domain d
         yrdpred<-mudpred+vd+ed            # Compute vector of out-if-sample responses
         ydnew<-c(ysd,yrdpred)             # Merge non-sample and sample values (full population or census) for domain d

         # Compute domain indicator using population values
         Ednew <- bxcx(ydnew,lambda=lambda,type=transform,InverseQ=TRUE)-m
         ebestimator[i,ell] <- indicator(Ednew)

      }  # End of Monte Carlo simulations for EB method

      ebestimator.EB[i]<-mean(ebestimator[i,])  # EB predictors of the measure (averages over the MC Monte Carlo simulations)

   } # End of cycle for target domain

   if (length_noselectdomXnonsample>0)
   {
      I <- length(selectdom_ini)
      ebestimator.EBaux <- ebestimator.EB
      ndaux <- nd
      ebestimator.EB <- rep(NA,I) 
      nd <- rep(0,I)
      for (i in 1:I)
      {
         pos <- which(selectdom==selectdom_ini[i])
         if (length(pos)!=0)
         {
            ebestimator.EB[i] <- ebestimator.EBaux[pos]
            nd[i] <- ndaux[pos]
         }
         else
            nd[i] <- sum(dom==d)        
      }
      selectdom <- selectdom_ini
   }   
   
   result$eb <- data.frame(domain=selectdom,eb=ebestimator.EB,sampsize=nd)
   result$fit <- Resultsfit

   return (result)

} #End of function 

