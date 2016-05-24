get.range.rLNLN <-
function(method,prev.result,dataset,n,TY,fix.mu,fixed.mu) {
# This function determines a parameter range in which it seems meaningful to search for
# a maximum of the log likelihood function, based on previous evaluations of the 
# log likelihood. The return value is a two-columnmatrix, where each row correspond to one 
# parameter; the first columns contains the lower bounds, the second column the upper bounds.
#
# The method can be
# - "none": no prior information used, return just NULL
# - "best": lower and upper bounds are the currently best estimates for each parameter
# - "1se": best parameter +/- 1 times se, where se is the empirical standard deviation of all entries
#   in prev.result
# - "2se": like "1se", just with +/- 2 times se
# - "top20": min and max of the top20 results
# - "quant": empirical 20% and 80% quantiles of all entries in prev.result
# - "mlci": approximate marginal 95% likelihood confidence intervals
      
   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)        
   stochprof.results <- NULL
   calculate.ci <- NULL        
   rm(stochprof.results)
   rm(calculate.ci)


   # trivial cases
   if (((method=="none") || (nrow(prev.result)==0)) || (is.null(prev.result))) {
      return(NULL)
   }
   
   # number of genes
   m <- (ncol(prev.result)+1)/TY - 2
   
   # results so far
   res <- stochprof.results(prev.result=prev.result,TY=TY,show.plots=F)
   if (is.null(res)) {
      return(NULL)
   }

   # best result so far   
   best <- res[1,-ncol(res),drop=F]
      
   ###################################
   # parameter ranges to be searched #
   ###################################
      
   if (method=="best") {
      # very best result; lower bound = upper bound
      ranges <- cbind(t(best),t(best))
   }
   else if (substr(method,2,3)=="se") {
      this.se <- apply(X=res[,-ncol(res)],MARGIN=2,FUN=sd)
      a <- as.double(substr(method,1,1))
      ranges <- cbind(t(best)-a*this.se,t(best)+a*this.se)
   }
   else if (method=="top20") {
      toptargets <- res[1:min(20,nrow(res)),1:(ncol(res)-1),drop=F]
      lower <- apply(X=toptargets,MARGIN=2,FUN=min)
      upper <- apply(X=toptargets,MARGIN=2,FUN=max)
      ranges <- cbind(lower,upper)
   }
   else if (method=="quant") {
      lower <- rep(NA,ncol(res)-1)
      upper <- lower
      for (i in 1:length(lower)) {
         this.par <- res[,i]
         this.par <- this.par[order(this.par)]
         lower.pos <- min(length(this.par),round(1,0.2 * length(this.par)+1))
         upper.pos <- max(1,round(0.8 * length(this.par)-1))
         lower[i] <- this.par[lower.pos]
         upper[i] <- this.par[upper.pos]   
      }
      ranges <- cbind(lower,upper)
   }
   else if (method=="mlci") {
      ranges <- calculate.ci(alpha=0.05,parameter=best,dataset=dataset,n=n,TY=TY,fix.mu=fix.mu,fixed.mu=fixed.mu)
      return(ranges)
   }
   
   ########################################
   # make sure the ranges above are valid #
   ########################################
      
   if (method!="best") {
      # p   
      if (TY>1) {
         ranges[1:(TY-1),1] <- pmax(0,ranges[1:(TY-1),1])
         ranges[1:(TY-1),2] <- pmin(1,ranges[1:(TY-1),2])
      }
      # sigma
      ranges[nrow(ranges)-((TY-1):0),1] <- pmax(ranges[nrow(ranges)-((TY-1):0),1],0.01)
      ranges[nrow(ranges)-((TY-1):0),2] <- pmax(ranges[nrow(ranges)-((TY-1):0),2],ranges[nrow(ranges)-((TY-1):0),1]+0.05)
   }
   
   return(ranges)
}
