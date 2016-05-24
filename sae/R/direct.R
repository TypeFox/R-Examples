direct <-
function(y,dom,sweight,domsize, data, replace=FALSE) {
    
   result <- data.frame(Domain=0,SampSize=0,Direct=0,SD=0,CV=0)

   missingsweight <- missing(sweight)
   missingdomsize <- missing(domsize)

  # direct estimator case
  # type 1: sampling without replacement
  # type 2: sampling without replacement under simple random sampling (SRS)
  # type 3: sampling with replacement
  # type 4: sampling with replacement under SRS

   if (replace==FALSE)
   { 
      if (!missingdomsize) 
      {
         if (!missingsweight)
            type <- 1
         else
            type <- 2
      } else
         stop("domsize is required when replace=FALSE.")
   } else if (replace==TRUE)
   { 
      if (!missingsweight)
      { 
         if (!missingdomsize)
            type <- 3
         else
            stop("domsize is required when replace=TRUE and sweight is used.")
      } else
         type <- 4
   } else
         stop("replace=",replace," must be TRUE or FALSE.")

    
#   classdata <- class(data)
#   if (classdata=="data.frame")
   if (!missing(data))
   {
      y       <- data[,deparse(substitute(y))]
      dom     <- data[,deparse(substitute(dom))]
      if (!missingsweight)
         sweight <- data[,deparse(substitute(sweight))]
   }

   if(!missingdomsize)
      if (any(is.na(domsize))) 
       stop(" domsize with NA values.")

   A<-length(y)
   B<-length(dom)
   if (!missingsweight)
   {
      C<-length(sweight)
      if((A!=B ) | (A!=C))
         stop(" y [",A,"], dom [",B,"] and sweight [",C,"] must be the same length.")
   } else
      if(A!=B)
         stop(" y [",A,"] and dom [",B,"] must be the same length.")


   # Delete rows with NA values
   rowNA <- c(which(is.na(y)),which(is.na(dom)))
   if (!missingsweight)  
      rowNA <- c(rowNA,which(is.na(sweight)))

   if (length(rowNA)>0)
   {
      y       <- y[-rowNA]
      dom     <- dom[-rowNA]
      if (!missingsweight)  
         sweight <- sweight[-rowNA]
   }

   # Test the domains are the same
   did     <- unique(dom)    # unique identifiers of domains
   Dsample <- length(did)    # number of domains in sample
   if (!missingdomsize)
   { 
      for (d in 1:Dsample)
      {
         ntimesdomi <- sum(did[d]==domsize[,1])
         if (ntimesdomi!=1)
            stop("Some sample domain indicators (dom) are not defined in population domain indicators.")
      }
   }
  

   # Calculate HT direct estimator for sampled domains   
   nds      <-rep(0,Dsample)   # domain sample sizes
   dirds    <-rep(0,Dsample)   # domain direct estimators 
   vardirds <-rep(0,Dsample)   # variances of direct estimators

   for (d in 1:Dsample) 
   {
      yd       <- y[dom==did[d]]
      nds[d]   <- length(yd)

      if (type==1)
      {
         sweightd <- sweight[dom==did[d]]        
         domsized <- domsize[(domsize[,1]==did[d]),2]

         dirds[d] <- sum(yd*sweightd)/domsized  
      
         # Approximated unbiased estimator of variance of HT direct estimator
         vardirds[d]<-sum(sweightd*(sweightd-1)*(yd^2))/(domsized^2)

      } else
      if (type==2)
      {
         domsized <- domsize[(domsize[,1]==did[d]),2]
         fd       <- nds[d]/domsized
         Sd2      <- var(yd)

         dirds[d]    <- sum(yd)/nds[d]  
         vardirds[d] <- (1-fd)*Sd2/nds[d]

      } else
      if (type==3)
      {
         sweightd <- sweight[dom==did[d]]
         domsized <- domsize[(domsize[,1]==did[d]),2]
         fd       <- nds[d]/domsized

         dirds[d]   <- sum(yd*sweightd)/domsized  
         vardirds[d]<- sum((fd*sweightd*yd-dirds[d])^2)/nds[d]

      } else  #type=4
      {
         Sd2   <- var(yd)

         dirds[d]    <- sum(yd)/nds[d]
         vardirds[d] <- Sd2/nds[d]
      } 
   }
  
   # direct estimator for non sampled domains
   if (!missingdomsize)
   {
      D <- nrow(domsize)
      if (D>Dsample)  
      {
         missing  <- rep(NA,D-Dsample)
         nd       <- c(nds,rep(0,D-Dsample)) # Domain direct estimators
         dird     <- c(dirds,missing)        # Domain direct estimators 
         vardird  <- c(vardirds,missing)     # Variances of direct estimators   
      
         domns    <- rep(0,D-Dsample)
         j<-1
         for (i in 1:D)
         {
            domi <- domsize[i,1]
            if (sum(domi==did)==0)
            {
               domns[j] <- domi
               j <- j+1
            }
         }
         did <- c(did,domns)
      } else
      {
         dird    <- dirds
         vardird <- vardirds
         nd      <- nds
      }
   } else
   {
      dird    <- dirds
      vardird <- vardirds
      nd      <- nds
   }

    
   # Percent coeficient of variation of direct estimator
   cvdird<-100*sqrt(vardird)/dird
  
   result   <- data.frame(Domain=did,SampSize=nd,Direct=dird,SD=sqrt(vardird),CV=cvdird)
   roworder <- order(result[,1])

   result   <- result[roworder,]
   return(result)	
}

