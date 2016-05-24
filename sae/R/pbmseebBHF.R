pbmseebBHF <-
function(formula,dom,selectdom,Xnonsample,B=100,MC=100,data,transform="BoxCox",lambda=0,constant=0,indicator)
{
   result <- list(est=NA, mse=NA)   

   if (missing(indicator))
      stop("Argument indicator must be included: indicator=function")

   if (transform!="BoxCox" & transform!="power")
      stop("Argument transform= \"",transform, "\" must be \"BoxCox\" or \"power\".")

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
      Xs  <- model.matrix(formula)        # Intercept is added if formula does not include -1
   }
   if (is.factor(dom))
      dom <- as.vector(dom)

   if (nrow(formuladataall)!=length(dom))
     stop("   Arguments formula [rows=",nrow(formuladataall),"] and dom [rows=",length(dom),"] must be the same length.\n")

   welfare <- formuladata[,1]            
   intercept <- (attributes(Xs)$assign[1]==0) 
   p    <- dim(Xs)[2]           # Number of auxiliary variables (included intercept if it is used)


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
      welfare  <- welfare[-rowNAdom]
      Xs  <- Xs[-rowNAdom,]
      dom <- dom[-rowNAdom]
   }

   # Domains for which EB estimators are required called (target domains)
   if (missing(selectdom))
      selectdom <- unique(dom)
   else
   {
      if (is.factor(selectdom))
         selectdom <- as.vector(selectdom)
      selectdom <- unique(selectdom)
   }

   if (intercept==TRUE)
      Xnonsample <- cbind(Xnonsample[,1],1,Xnonsample[,-1])

   ## Fit the model
   result$est <- ebBHF(welfare~Xs-1,dom=dom,selectdom=selectdom,Xnonsample=Xnonsample,MC=MC,    
                       transform=transform,lambda=lambda,constant=constant,indicator=indicator)  

   betaest <- matrix(result$est$fit$fixed,nrow=p,ncol=1) # Vector of model coefficients (size p) 
   sigmae2est <- result$est$fit$errorvar          # Estimated error variance
   sigmau2est <- result$est$fit$refvar            # VarCorr(fit2) is the estimated cov. matrix of the model random components


   I    <- length(selectdom)
   n    <-length(welfare)  # Total number of sample observations
   udom <- unique(dom)
   D    <-length(udom)     # Total number of domains with sample data

   ## Test that domain codes in selectdom are defined in Xnonsample
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

   # Sample sizes and out-of-sample sizes of target domains
   nd <- rep(0,D) 
   rd <- rep(0,I)
   musd.B <- murd.B <- list()

   for (d in 1:D)
   {
      domd    <- udom[d]
      rowsd   <- (dom==domd)
      nd[d]   <- sum(rowsd)

      posdomd  <- which(udom==domd) 
      musd.B[[posdomd]] <- Xs[rowsd,]%*%betaest       # Take sample elements
   }

   for (i in 1:I)
   {
      d    <- selectdom[i]                   
      Xaux <- as.matrix(Xnonsample[(Xnonsample[,1]==d),-1])                     # delete first column with domain codes
      murd.B[[i]] <- Xaux%*%betaest   # Vector of out-of-sample marginal means
      rd[i]  <- dim(murd.B[[i]])[1]
   }

  # Transform welfare after adding a constant to make it positive
   #m <- max(c(0,(-1)*min(welfare)+1))  
   m <- constant           
   
   MSEindicatorsum.B<-rep(0,I)  # Initialize vectors with the bootstrap MSEs of EB estimators
   MSEindicatorEB.B <-rep(0,I)
#   trueindicator.B  <-rep(0,I)  # Initialize vectors with true values of parameters
#   indicator.EB.B   <-rep(0,I)  # Initialize vectors with estimated poverty measures

   cat("\nBootstrap procedure with B =",B,"iterations starts.\n")

   #set.seed(123)  # prueba

   sqrtsigmae2est <- sqrt(sigmae2est)
   sqrtsigmau2est <- sqrt(sigmau2est)

   warn_ini <- getOption("warn") # Do not show eblupBHF warning messages  
   options(warn = -1)

   for (b in 1:B)        ### Start of bootstrap cycle
   {   
      cat("b =",b,"\n")  

      # Generate a bootstrap population by generating sample and 
      # out-of-sample elements and calculate values of the parameters.

      ys.B<-rep(0,n)   # Bootstrap sample vector
      ud.B<-rep(0,D)   # Bootstrap random effects
      
      for (d in 1:D)         # Generate sample elements (for all domains)
      {
         domd    <- udom[d]                        
         posdomd <- which(udom==domd)
         rowsd   <- (dom==domd)    

         # Generate sample values of y from the fitted model
         esd.B         <- rnorm(nd[posdomd],0,sqrtsigmae2est)  
         ud.B[posdomd] <- rnorm(1,0,sqrtsigmau2est)
         ys.B[rowsd]   <- musd.B[[posdomd]]+ud.B[posdomd]+esd.B
      }

      # Compute EB predictors for each bootstrap sample
      Es.B <- bxcx(ys.B,lambda=lambda,type=transform,InverseQ=TRUE)-m  
      resultsEB<-ebBHF(Es.B~Xs-1,dom=dom,selectdom=selectdom,Xnonsample=Xnonsample,MC=MC,
                         transform=transform,lambda=lambda,constant=constant,indicator=indicator)
  
      # Generate non-sample values only for selected domains
      for (i in 1:I)
      {
         d      <- selectdom[i]
         posd   <- which(udom==d)
         erd.B  <- rnorm(rd[i],0,sqrt(sigmae2est))   # Vector of out-of-sample random errors in target domain i

         if (length(posd)==0){ ###DUDA
            # Vector of out-of-sample responses
            yrd.B <- murd.B[[i]]+rnorm(1,0,sqrtsigmau2est)+erd.B #DUDA

            Erd.B <- bxcx(yrd.B,lambda=lambda,type=transform,InverseQ=TRUE)-m  
            Ed.B  <- Erd.B
         } else
         {
            # Vector of out-of-sample responses
            yrd.B <- murd.B[[i]]+ud.B[posd]+erd.B 
            Erd.B <- bxcx(yrd.B,lambda=lambda,type=transform,InverseQ=TRUE)-m  
 
           # Merge sample and non-sample elements for target domain d
            ysd.B <- ys.B[dom==d]
            Esd.B <- bxcx(ysd.B,lambda=lambda,type=transform,InverseQ=TRUE)-m  
            Ed.B  <- c(Erd.B,Esd.B)
         }
        
#         trueindicator.B[i] <- indicator(Ed.B)             # Calculate true parameter for target domain
#         indicator.EB.B[i]  <-resultsEB$EB$ebestimator[i]  # EB predictors of the measure
#         MSEindicatorsum.B[i] <- MSEindicatorsum.B[i] +(indicator.EB.B[i]-trueindicator.B[i])^2   # Cumulated squared errors for Bootstrap MSE
         MSEindicatorsum.B[i] <- MSEindicatorsum.B[i] +(resultsEB$eb$eb[i]-indicator(Ed.B))^2   # Cumulated squared errors for Bootstrap MSE                                                                                                

       }
   }  # End of the bootstrap cycle

   options(warn = warn_ini)

   # Bootstrap MSS of EB estimators
   MSEindicatorEB.B <- MSEindicatorsum.B/B

   if (length_noselectdomXnonsample>0)
   {
      I     <- length(selectdom_ini)
      pbmse <- rep(NA,I)
      for (i in 1:I)
      {
         pos <- which(selectdom==selectdom_ini[i])
         if (length(pos)!=0)
            pbmse[i] <- MSEindicatorEB.B[pos]
      }
      selectdom <- selectdom_ini
   } else
      pbmse <- MSEindicatorEB.B

   result$mse <-data.frame(domain=selectdom,mse=pbmse)
   return (result)

} # End of function PBMSE.EB
