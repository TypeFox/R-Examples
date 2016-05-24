stochprof.results.LNLN <-
function(prev.result,TY,show.plots=T,plot.title="",pdf.file,fix.mu=F) {
# Evaluates the set of results that are contained in prev.result. That means, it removes
# entries where the target function is equal to infinity, it removes double entries,
# it removes unlikely parameter combinations (if there are too many) etc. When show.plots==T, 
# the results are graphically displayed.
#
# Parameters:
# - prev.result contains parameter combinations and the respective value of the target 
#   function. It is typically the output of "stochprof.search".
# - TY is the number of types of cells assumed in the model.
# - If show.plots==T, the results are plotted. In particular, one plot is produced for each
#   parameter, with the value of the parameter plotted against the value of the target 
#   function. This is not exactly the profile log-likelihood function because there is
#   no conditioning on the other parameters being equal to the ML estimate. If the estimation
#   procedure has converged, however, one can recognize the shape of the profile log-likelihood
#   from these plots. A red bar indicates the position of the maximum likelihood estimator.
# - If show.plots==T, plot.title (default="") is the title of each plot.
# - If pdf.file is not missing, the plots will be written into this pdf.file. The file has
#   to include the entire path!
# - fix.mu (default=F) indicates whether the log-mean has been kept fixed. In that case, no
#   plots will be produced for these parameters.

   # are there actually any results in prev.result?
   if (is.null(prev.result) || (nrow(prev.result)==0)) {
      print("stochprof.results: File contains only the header.")
      return(NULL)
   }

   # "results" will be the return variable
   results <- prev.result
   these.names <- colnames(results)      

   # number of genes
   m <- (ncol(results)-1)/TY - 1
      
   ###############
   # sorting out #
   ###############        
        
   # remove those rows with infinite target values
   inf.indices <- which(results[,"target"]>=10^7)
   if (length(inf.indices)>0) {
      results <- results[-inf.indices,,drop=F]
   }  
   if (nrow(results)==0) {
      print("stochprof.results: dataset contains only infinite target values.")
      return(NULL)
   }
            
   # remove those rows with very small or very large mu1, mu2
   mu.min <- apply(X=results[,TY:((m+1)*TY-1),drop=F],FUN=min,MARGIN=1)
   mu.max <- apply(X=results[,TY:((m+1)*TY-1),drop=F],FUN=max,MARGIN=1)      
   smallmu.indices <- which(mu.min<(-50))
   largemu.indices <- which(mu.max>10)
   mu.indices <- union(smallmu.indices,largemu.indices)
   if ((length(mu.indices)>0) && (length(mu.indices)<nrow(results))) {
      results <- results[-mu.indices,,drop=F]
   }
             
   #########################
   # remove double entries #
   #########################

   # round all entries to these numbers of decimals            
   nod.p <- 4
   nod.all <- 3
   nod.target <- 4
      
   if (TY>1) {
      p.indices <- 1:(TY-1)
   }
   else {
      p.indices <- NULL
   }
   target.indices <- (m+1)*TY+1
            
   results[,p.indices] <- round(results[,p.indices],nod.p)   
   results[,-c(p.indices,target.indices)] <- round(results[,-c(p.indices,target.indices)],nod.all)
   results[,target.indices] <- round(results[,target.indices],nod.target)  

   # remove doubles
   results <- unique(results,MARGIN=1)      

   ########    
   # sort #
   ########
   
   # small target (i.e. small negative log likelihood) is the best; 
   # choose ascending order
   results <- results[order(results[,"target"]),,drop=F]
                  
      
   ###################################
   # further reduce number of points #
   ###################################
   
   indices <- which(results[,"target"]<=log(5) + (results[,"target"])[1])

   minlength <- 50      
   maxlength <- 5000              

   # if there are less than "minlength" entries in "indices", keep this minimum number
   if (length(indices)<minlength) {
      indices <- 1:min(minlength,nrow(results))
   }
   # if there are more than "maxlength" entires in "indices", keep only this maximum number
   else if (length(indices)>maxlength) {
      indices <- 1:min(maxlength,nrow(results))
   }
                
   results <- results[indices,,drop=F]

   
   #########
   # plots #
   #########
   if (show.plots) {
      # mark the estimate with smallest target
      best <- results[1,,drop=F]
     
      # which variables to plot

      # by default: all variables except target
      index.set <- 1:(ncol(results)-1)
      # but if fix.mu==T, remove those indices
      if (fix.mu) {        
         index.set <- index.set[-(TY:((m+1)*TY-1))]
      }
      
      if (missing(pdf.file)) {
         par(ask=T)
      }
      else {
         pdf(pdf.file)
      }

      for (i in index.set) {
         plot(results[,i],-results[,"target"],xlab=these.names[i],ylab="log likelihood",main=plot.title,cex.axis=1.5,cex.lab=1.5,pch=1,lwd=3)
         abline(v=best[i],lwd=2,col="red")
      }

      if (missing(pdf.file)) {      
         par(ask=F)
      }
      else {
         dev.off()
      }
   }
   
   ################
   # return value #
   ################
   
   return(results)
}
