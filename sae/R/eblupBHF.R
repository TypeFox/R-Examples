eblupBHF <-
function(formula,dom,selectdom,meanxpop,popnsize,method="REML",data)
{
   result <- list(eblup=NA, fit=NA) 

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


   ys <- formuladata[,1]            
   intercept <- (attributes(Xs)$assign[1]==0) 
   p<-dim(Xs)[2]                   # Number of auxiliary variables (including the constant)


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
      ys  <- ys[-rowNAdom]
      Xs  <- Xs[-rowNAdom,]
      dom <- dom[-rowNAdom]
   }
  
# Number of domains for which EBLUPs are required called (target domains)
   if (missing(selectdom))
      selectdom <- unique(dom)
   else
   {
      if (is.factor(selectdom))
         selectdom <- as.vector(selectdom)
      selectdom <- unique(selectdom)
   }
   I <- length(selectdom)
   n <- length(ys)                # Total number of sample observations

# Save meanxpop and popnsize domains in the same order as selected domains
# If a selected domain is not in meanxpop or popnsize domains, the value is 0
# and a warning is printed

   meanxpopsel           <- data.frame(matrix(0,nrow=I,ncol=ncol(meanxpop)))
   colnames(meanxpopsel) <- colnames(meanxpop)
   popnsizesel           <- data.frame(matrix(0,nrow=I,ncol=ncol(popnsize)))
   colnames(popnsizesel) <- colnames(popnsize)
   
   noselectdom_meanxpop <- NULL
   noselectdom_popnsize <- NULL
   selectdomfinded      <- NULL
   for (i in 1:I)
   {
      idom1 <- which(selectdom[i]==meanxpop[,1])
      idom2 <- which(selectdom[i]==popnsize[,1])
      if (length(idom1)!=0 && length(idom2)!=0)  #if selectdom finded
      {
         meanxpopsel[i,] <- meanxpop[idom1,]
         popnsizesel[i,] <- popnsize[idom2,]
         selectdomfinded <- c(selectdomfinded,selectdom[i])
      } else
      {
         if (length(idom1)==0)
            noselectdom_meanxpop <- c(noselectdom_meanxpop,selectdom[i])
         if (length(idom2)==0)
            noselectdom_popnsize <- c(noselectdom_popnsize,selectdom[i])
      }
   }
   ## warning message
   length1 <- length(noselectdom_meanxpop)
   length2 <- length(noselectdom_popnsize)
   if (length1>0 || length2>0 )
   {
      text <- "The following domain codes (selectdom) are not defined in population domain codes."
      text1 <- text2 <- ""
      if (length1>0)
      {
         noselectdom_meanxpop <- sort(noselectdom_meanxpop)
         for (i in 1:length1)
            text1 <- paste(text1,noselectdom_meanxpop[i])
         text <- paste(text, "\n    - meanxpop[,1]:", text1,sep="")
      }
      if (length2>0)
      {
         noselectdom_popnsize <- sort(noselectdom_popnsize)
         for (i in 1:length2)
            text2 <- paste(text2,noselectdom_popnsize[i])
         text <- paste(text, "\n    - popnsize[,1]:", text2,sep="")
      }   
      warning(text)
   }

   meanxpop <- as.matrix(meanxpopsel[,-1])
   popnsize <- as.matrix(popnsizesel[,-1])

   if ((intercept==TRUE) && p==ncol(meanxpop)+1)             
      meanxpop <- cbind(1,meanxpop)

# Fit the nested-error model to sample data by REML/ML 
# method using function lme from library nlme. 
   fit.EB<-try(lme(ys~-1+Xs,random=~1|as.factor(dom),method=method))
   if (class(fit.EB)=="try-error")   # Aqui en el caso de que se llame desde pbmseBHF
                                     # habria que generar otra muestra
   {
      #options(warn = 0)
      #stop("lme function within eblupBHF.")
      warning("lme function within eblupBHF.")
      return (result)
   }

# Create a list object containing different results from the model fit.
   Resultsfit<-list(summary=summary(fit.EB), fixed=fixed.effects(fit.EB),
                    random=random.effects(fit.EB), errorvar=fit.EB$sigma^2,
                    refvar=as.numeric(VarCorr(fit.EB)[1,1]), loglike=fit.EB$logLik,
                    residuals=fit.EB$residuals[1:n])

# Save some of the results of the fitting method in variables
   betaest   <-matrix(fixed.effects(fit.EB),nrow=p,ncol=1)  # Vector of model coefficients (size p)
   upred     <-random.effects(fit.EB)                   # Predicted random effects: Watch out! It is not a vector, it is a matrix with 1 column

   eblup <- rep(0,I)
   meanXsd <- matrix(0,nrow=1,ncol=p)
   text <- NULL
   SampSizeselectdom <- rep(0,I)
   for (i in 1:I)          # Cycle for target domains 
   {        
     d  <- selectdom[i]    # target domain code
     if (length(which(d==selectdomfinded))==0)
        eblup[i]<-NA
     else
     {
        rowsd <- (dom==d)        # Get sample values for target domain
        if (any(rowsd)==TRUE)
        {
           Xsd   <- matrix(Xs[rowsd,],ncol=p)
           fd    <- sum(rowsd)/popnsize[i]
           for (k in 1:p) 
              meanXsd[1,k]<-mean(Xsd[,k])
           eblup[i]<-fd*mean(ys[rowsd])+(meanxpop[i,]-fd*meanXsd)%*%betaest+(1-fd)*upred[as.character(d),1]
           SampSizeselectdom[i] <- sum(rowsd)
        } else
        {
           eblup[i]<- meanxpop[i,]%*%betaest    ## PREGUNTAR A ISABEL
           text    <- paste(text,d)
        } 
      }
   }  

   if (is.null(text)==FALSE)
      warning("The following selected domains (selectdom) have zero sample size. \n  The EBLUPs of these domains are the synthetic regression estimators.\n   Domains: ",text)  

   result$eblup <- data.frame(domain=selectdom,eblup=eblup,sampsize=SampSizeselectdom)
   result$fit   <- Resultsfit 
   return (result)

} #End of function

