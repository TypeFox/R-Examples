pssynt <-
function(y,sweight,ps,domsizebyps,data) 
{
   result <- data.frame(Domain=0,PsSynthetic=0)

   if (!missing(data))
   {
      y <- data[,deparse(substitute(y))]
      sweight <- data[,deparse(substitute(sweight))]
      ps <- data[,deparse(substitute(ps))]
   }

   if (any(is.na(domsizebyps))) 
       stop(" domsizebyps with NA values.")

   A<-length(y)
   B<-length(sweight)
   C<-length(ps)
   if((A!=B ) | (A!=C))
         stop(" y [",A,"], sweight [",B,"] and ps [",C,"] must be the same length.")

   # Delete rows with NA values
   rowNA <- c(which(is.na(y)),which(is.na(sweight)),which(is.na(ps)))
   if (length(rowNA)>0)
   {
      y       <- y[-rowNA]
      sweight <- sweight[-rowNA]
      ps      <- ps[-rowNA]
   }


   # Test the ps sample identifiers are in ps population identifiers  
   psidsample <- unique(ps)                            # Unique poststrata identifiers in sample
#   psidpop    <- as.integer(colnames(domsizebyps)[-1]) # Poststrata identifiers in population
   psidpop    <- (colnames(domsizebyps)[-1]) # Poststrata identifiers in population

   Jsample    <- length(psidsample)      # Number of poststrata in sample
   idpssamplenofinded <- NULL
   for (j in 1:Jsample)
   {
      find <- which(psidsample[j]==psidpop)
      if (length(find)==0)
         idpssamplenofinded <- c(idpssamplenofinded,psidsample[j])
   }
   if (length(idpssamplenofinded)>0)
   {
      text1 <- text2 <- ""

      idpssort <- sort(idpssamplenofinded)
      for (i in 1:length(idpssamplenofinded))
         text1 <- paste(text1,idpssort[i])

      idpssort <- sort(psidpop)
      for (i in 1:length(idpssort))
         text2 <- paste(text2,idpssort[i])
         
      warning("Some post-strata indicators (ps) are not defined in population post-strata indicators (domsizebyps[,1]).",
              "\n    - ps not defined:",text1,
              "\n    - population ps :",text2
              )
   }
 
   Ndj     <- as.matrix(domsizebyps[,-1]) # popn. sizes of each intersection of domain and post-strata
   Nj      <- colSums(Ndj)                # popn sizes of each post-strata j
   domsize <- rowSums(Ndj)                # popn sizes of each domain
  
   J    <- ncol(Ndj)
   dirj <- matrix(0,nrow=1,ncol=J)        # Horvitz-Thompson direct estimators of the post-strata means
   colnames(dirj) <- colnames(Ndj) <- psidpop

   for (j in 1:Jsample)
   {
      psid       <- psidsample[j]
      labelpsid  <- as.character(psid) 
      if (length(which(labelpsid==idpssamplenofinded))==0)
      {
         yj         <- y[ps==psid]
         sweightj   <- sweight[ps==psid]
         dirj[1,labelpsid] <- sum(yj*sweightj)/Nj[labelpsid]
      }
   }

  # Post-stratified estimators for each domain. 
   pssynt   <- as.vector(Ndj%*%t(dirj))/domsize 
  
   result   <- data.frame(Domain=domsizebyps[,1],PsSynthetic=pssynt)
   roworder <- order(result[,1])

   result   <- result[roworder,]
   return(result)	
}
