ssd <-
function(dom,sweight,domsize,direct,synthetic,delta=1,data) {

   result <- data.frame(Domain=0,ssd=0,CompWeight=0)
   
   if (!missing(data))
   {
      dom <- data[,deparse(substitute(dom))]
      sweight <- data[,deparse(substitute(sweight))]
   }

   if (any(is.na(domsize))) 
       stop(" domsize with NA values.")

   # dom, sweight the same lenght
   A <- length(dom)
   B <- length(sweight)
   if(A!=B)
      stop(" dom [",A,"] and sweight [",B,"] must be the same length.")

   # Check direct, synthetic and domsize have the same length (D) 
   nrowdirect    <-nrow(direct)
   nrowsynthetic <-nrow(synthetic)
   nrowdomsize   <-nrow(domsize)
   if((nrowdirect!=nrowsynthetic ) | (nrowdirect!=nrowdomsize))
      stop("direct [",nrowdirect,"] synthetic [",nrowsynthetic,"] and domsize [",nrowdomsize,"] must be the same length")

   # Order direct, synthetic and domsize by domains in ascending order
   roworderdirect    <- order(direct[,1])
   rowordersynthetic <- order(synthetic[,1])
   roworderdomsize   <- order(domsize[,1])

   direct    <- direct[roworderdirect,]
   synthetic <- synthetic[rowordersynthetic,]
   domsize   <- domsize[roworderdomsize,]

   # Test the domains are the same 
   equaldomains <- sum(direct[,1]==synthetic[,1])
   if (equaldomains!=nrowdirect)
      stop("direct and synthetic have some different domain indicators.")

   equaldomains <- sum(direct[,1]==domsize[,1])
   if (equaldomains!=nrowdirect)
      stop("direct and domsize have some different domain indicators.")

   iddom   <- domsize[,1]      # domain identifiers in ascending order   
   did     <- unique(dom)      # unique identifiers of sampled domains

   Dsample  <- length(did)
   for (d in 1:Dsample)
   {
      ntimes <- sum(did[d]==iddom)
      if (ntimes!=1)
         stop("Some sample domains indicators (dom) are not defined in population domain indicators.")
   }


   D     <- nrowdomsize # Number of domains in population
   Ndhat <- rep(0,D)    # Estimated domain population sizes
   phid  <- rep(1,D)    # Weights attached to direct estimators in the composition
   for (d in 1:D) 
   {
        Ndhat[d]     <- sum(sweight[dom==iddom[d]])
        deltadomsize <- delta*domsize[d,2]
        if (Ndhat[d]<deltadomsize)
           phid[d] <- Ndhat[d]/deltadomsize
   }

   # Calculate SSD estimator
   ssd<-phid*direct[,2]+(1-phid)*synthetic[,2]
    
   return(data.frame(Domain=iddom,ssd=ssd,CompWeight=phid))
}

