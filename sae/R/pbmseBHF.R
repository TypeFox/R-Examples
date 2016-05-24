pbmseBHF <-
function(formula,dom,selectdom,meanxpop,popnsize,B=200,method="REML",data)
{
   result <- list(est=NA, mse=NA)

   if (!missing(data))
   {
      formuladataall <- model.frame(formula,na.action = na.pass,data)
      formuladata <- model.frame(formula,na.action = na.omit,data)
      Xs  <- model.matrix(formula,data)        # Intercept is added if formula does not include -1
      dom <- data[,deparse(substitute(dom))]
   } else
   {
      formuladataall <- model.frame(formula,na.action = na.pass)
      formuladata <- model.frame(formula,na.action = na.omit)
      Xs  <- model.matrix(formula)             # Intercept is added if formula does not include -1
   }
   if (is.factor(dom))
      dom <- as.vector(dom)

   if (nrow(formuladataall)!=length(dom))
     stop("   Arguments formula [rows=",nrow(formuladataall),"] and dom [rows=",length(dom),"] must be the same length.\n")

   ys <- formuladata[,1]            
   intercept <- (attributes(Xs)$assign[1]==0) 
   p<-dim(Xs)[2]                           # Number of auxiliary variables (including the constant)


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


   if ((intercept==TRUE) && p==ncol(meanxpop))             
      meanxpopaux <- cbind(meanxpop[,1],1,meanxpop[,-1])
   else
      meanxpopaux <- meanxpop

# Fit the nested-error model to sample data by REML/ML method using function lme from library nlme. 
   result$est <- eblupBHF(ys~Xs-1,dom,selectdom,meanxpopaux,popnsize,method)
   if (is.data.frame(result$est$eblup) == FALSE)
   {
      warning("The fitting method does not converge.\n") 
      return (result);
   }


   betaest    <- matrix(result$est$fit$fixed,nrow=p,ncol=1) # Vector of model coefficients (size p) 
   sigmae2est <- result$est$fit$errorvar          # Estimated error variance
   sigmau2est <- result$est$fit$refvar            # VarCorr(fit2) is the estimated cov. matrix of the model random components


   I    <- length(selectdom)
   n    <- length(ys)           # Total number of sample observations


# Save meanxpop and popnsize domains in the same order as selected domains
# If a selected domain is not in meanxpop or popnsize domains, the value is 0
# and a warning is printed

   meanxpopsel           <- data.frame(matrix(0,nrow=I,ncol=ncol(meanxpop)))
   colnames(meanxpopsel) <- colnames(meanxpop)
   popnsizesel           <- data.frame(matrix(0,nrow=I,ncol=ncol(popnsize)))
   colnames(popnsizesel) <- colnames(popnsize)
   
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
      } 
   }

   meanxpop <- as.matrix(meanxpopsel[,-1])  # Delete domain codes
   popnsize <- as.matrix(popnsizesel[,-1])

   if ((intercept==TRUE) && p==ncol(meanxpop)+1)             
      meanxpop <- cbind(1,meanxpop)
##################

   udom   <- unique(dom)
   D      <- length(udom)         # Total number of domains in sample data
   nd     <- rep(0,D)
   SampSizeselectdom <- rep(0,I)
   musd.B <- mud.B <- list() #

   for (d in 1:D)
   {
      rowsd       <- (dom==udom[d])    
      musd.B[[d]] <- Xs[rowsd,]%*%betaest  # Take sample elements
      nd[d]       <- sum(rowsd)            # Sample sizes of sampled domains
   }

   for (i in 1:I)
   {
       mud.B[[i]]<-sum(meanxpop[i,]*betaest[,1]) #
#      murd.B[[i]] <- Xnonsample[[i]]%*%betaest
                                           # Sample sizes of selected domains
      d    <- selectdom[i]
      posd <- which(udom==d)
      if (length(posd)>0)
         SampSizeselectdom[i] <- nd[posd]
   }      

   Ni <- popnsize
   rd <- Ni-SampSizeselectdom

   meanxpop <- data.frame(selectdom,meanxpop)
   popnsize <- data.frame(selectdom,popnsize)
   
   MSE.B    <- truemean.B <- rep(0,I)      # Initialize Vectors with the bootstrap MSEs of EBLUPs

   warn_ini <- getOption("warn") # Do not show eblupBHF warning messages  
   options(warn = -1)

   cat("\nBootstrap procedure with B =",B,"iterations starts.\n")
   b <- 1
   while (b<=B)   ### Bootstrap cycle starts
   {   
      ys.B      <- rep(0,n)      # Bootstrap sample vector
      ud.B      <- rep(0,D)      # Bootstrap random effects
      esdmean.B <- rep(0,D)
      for (d in 1:D)             # Generate sample elements (for sampled areas)
      {
         esd.B        <- rnorm(nd[d],0,sqrt(sigmae2est))  
         ud.B[d]      <- rnorm(1,0,sqrt(sigmau2est))
         rowsd        <- (dom==udom[d]) 
         ys.B[rowsd]  <- musd.B[[d]]+ud.B[d]+esd.B    # Generate y sample values from the fitted model
         esdmean.B[d] <- mean(esd.B) #
      }

      # Calculate bootstrap truemeans for selected domains
      for (i in 1:I)
      {
         erdmean.B <- rnorm(1,0,sqrt(sigmae2est/rd[i])) #
         posd      <- which(udom==selectdom[i])
         if (length(posd)!=0)
         {
            edmean.B     <-esdmean.B[posd]*nd[posd]/Ni[i]+erdmean.B*rd[i]/Ni[i] #
            truemean.B[i]<-mud.B[[i]]+ud.B[posd]+edmean.B #
         } else
            truemean.B[i]<-mud.B[[i]]+rnorm(1,0,sqrt(sigmau2est))+erdmean.B #
      } 

      # Compute EBLUPs for each bootstrap sample
      mean.EB <- eblupBHF(ys.B~Xs-1,dom=dom,selectdom=selectdom,meanxpop=meanxpop,popnsize=popnsize,method=method)  
      if (is.data.frame(mean.EB$eblup) == FALSE)
         next

      cat("b =",b,"\n")

      MSE.B   <- MSE.B+(mean.EB$eblup$eblup-truemean.B)^2 # Cummulated squared errors for Bootstrap MSE
      b <- b+1    
   }# End of bootstrap cycle 

   options(warn = warn_ini)

   MSEEB.B<-MSE.B/B     # Bootstrap MSE of EB estimators

# Creo que se podria quitar,... COMPROBAR
#   text <- NULL
#   for (i in 1:I)          # Cycle for target domains 
#   {        
#     d  <- selectdom[i]    # target domain code
#     if ((SampSizeselectdom[i]==0) && (length(which(d==selectdomfinded))!=0))
#        text <- paste(text,d)
#   }  
#   if (is.null(text)==FALSE)
#      warning("The following selected domains (selectdom) have zero sample size. \n  The EBLUPs of these domains are the synthetic regression estimators.\n   Domains: ",text)  

   result$mse <-data.frame(domain=selectdom,mse=MSEEB.B)
   return (result)

} # End of function 
