desc.wts <-function(data,w, sampw = sampw, 
                    vars=NULL,
                    treat.var,
                    tp,
                    na.action="level",
                    perm.test.iters=0,
                    verbose=TRUE,
                    alerts.stack,
                    estimand, multinom = FALSE, fillNAs = FALSE){
   if(is.null(vars)) vars <- names(data)[names(data)!=treat.var]
   ess.ctrl          <- (sum(w[data[,treat.var]==0])^2)/sum(w[data[,treat.var]==0]^2)
   ess.treat    <- (sum(w[data[,treat.var]==1])^2)/sum(w[data[,treat.var]==1]^2)
   
   sampW <- sampw
   vars1 <- vars
   treat.var1 <- treat.var


#multinom <- FALSE

   bal.tab   <- bal.stat(data=data,w.all=w, sampw = sampW,
                         vars=vars1,
                         treat.var=treat.var1,
                         na.action=na.action,
                         estimand=estimand, multinom = multinom, fillNAs = fillNAs)
   pval.maxks <- NA
   # compute permutation p-values for KS statistic
   if(perm.test.iters>0)
   {
      ess.t <- (sum(w[data[,treat.var]==1])^2)/
               sum(w[data[,treat.var]==1]^2)
      ess.c <- (sum(w[data[,treat.var]==0])^2)/
               sum(w[data[,treat.var]==0]^2)
      w.1   <- w

   ### revised 091710
      w.1[data[,treat.var]==0] <- ess.c*w[data[,treat.var]==0]/sum(w[data[,treat.var]==0])
      w.1[data[,treat.var]==1] <- ess.t*w[data[,treat.var]==1]/sum(w[data[,treat.var]==1])
      w.1 <- w.1/(ess.c+ess.t)

      pval.var   <- rep(0,nrow(bal.tab$results))
      names(pval.var) <- rownames(bal.tab$results)
      pval.maxks <- 0
      if(verbose)
      {
         cat("Permutation test progress:\n")
         progress <- round(perm.test.iters*c(0.01,(1:10)/10))
      }
      for(i.rep in 1:perm.test.iters)
      {
         if(verbose && (i.rep %in% progress))
         {
            cat(round(100*i.rep/perm.test.iters),"%\n",sep="")
         }
         i <- sample(1:nrow(data),ess.t+ess.c,replace=TRUE,prob=w.1)
         temp <- data[i,c(vars,treat.var)]
         temp[,treat.var] <- as.numeric((1:nrow(temp))<=ess.t)
         bal.temp   <- bal.stat(data=temp,
                                w.all=rep(1,length(i)),
                                sampw = rep(1,length(i)),
                                vars=vars1,
                                treat.var=treat.var1,
                                na.action=na.action,
                                get.means=FALSE,
                                estimand=estimand, multinom = multinom)
         # if a bootstrap sample is missing a level need to make sure
         #   these still align
         j <- match(rownames(bal.temp$results),names(pval.var))
         pval.var[j] <- pval.var[j]   + 
                        (bal.temp$results$ks >= bal.tab$results$ks[j])
         pval.maxks  <- pval.maxks + 
                        (max(bal.temp$results$ks) >= max(bal.tab$results$ks))
      }
      pval.var   <- as.numeric(pval.var/perm.test.iters)
      pval.maxks <- as.numeric(pval.maxks/perm.test.iters)
      
      # replace the analytic p-values computed in bal.tab
      bal.tab$results$ks.pval <- pval.var
   }
   
   check.err(cov.table=bal.tab$results, stage=tp, alerts.stack=alerts.stack, estimand=estimand, ess.ctrl=ess.ctrl, ess.treat=ess.treat)

   max.ks  <- max(bal.tab$results$ks)
   mean.ks  <- mean(bal.tab$results$ks)
   max.es <- with(bal.tab$results, max(abs(std.eff.sz[std.eff.sz<500]),
                 na.rm=TRUE))
   mean.es   <- with(bal.tab$results, mean(abs(std.eff.sz[std.eff.sz<500]),
                  na.rm=TRUE))
   return(list(ess.ctrl=ess.ctrl,
               ess.treat=ess.treat,
               n.treat=sum(data[,treat.var]==1),
               n.ctrl =sum(data[,treat.var]==0),
               max.es=max.es,
               mean.es=mean.es,
               max.ks=max.ks,
               max.ks.p=pval.maxks,
               mean.ks=mean.ks,
               bal.tab=bal.tab))
}

