stochprof.loop <-
function(model,dataset,n,TY,genenames=NULL,fix.mu=F,fixed.mu,par.range=NULL,prev.result=NULL,loops=5,until.convergence=T,print.output=T,show.plots=T,plot.title="",pdf.file,use.constraints=F,subgroups) {
# This function carries out maximum likelihood estimation for the parameters in the 
# stochastic profiling model. Because the log-likelihood function is multimodal, no
# straightforward use of gradient-based approaches for finding globally optimal parameter
# combinations is possible. To tackle this challenge, this function performs a two-step estimation 
# procedure: First, it computes the log-likelihood function at randomly drawn parameter
# combinations to identify high-likelihood regions in parameter space at computationally
# low cost. Then, it uses the Nelder-Mead algorithm to identify local maxima of the
# likelihood function. The starting values for this algorithm are randomly drawn from the
# high-likelihood regions determined in the first step. To further localize the global
# optimum, the function again performs grid searches of the parameter space, this time 
# around the optimum identified by the Nelder-Mead algorithm. This search creates another 
# space to identify high-likelihood regions, which are then used to seed another Nelder-Mead
# optimization.
#
# Parameters:
#   
# - model is the chosen model, either "LN-LN", "rLN-LN" or "EXP-LN"
# - dataset is a matrix which contains the cumulated expression data over all cells in a tissue sample.
#   Columns represent different genes, rows represent different tissue samples.
# - n is the number of cells taken from each tissue sample.
# - TY is the number of types of cells that is assumed in the stochastic model.
# - genenames (default=NULL) are the names of the genes in the dataset.
#   For genenames==NULL, the genes will simply be enumerated according to the column numbers in the dataset.
# - fix.mu (default=F) indicates whether the log-means are kept fixed in the estimation procedure or whether
#   they are to be estimated.
# - fixed.mu (no default, needs to be specified only when fix.mu==T) is a vector containing the values to which
#   the log-means should be fixed if fix.mu==T. The order of components is as follows:
#   (mu_type_1_gene_1, mu_type_1_gene_2, ..., mu_type_2_gene_1, mu_type__gene_2, ...)
# - par.range (default=NULL) is a range from which the parameter values should be randomly drawn if there is no 
#   other knowledge available. This is a matrix with the number of rows being equal to the number of model 
#   parameters. The first column contains the lower bound, the second column the upper bound. If par.range==NULL, 
#   some rather large range is defined.
# - prev.result (default=NULL) can contain results from former calls of this function.
# - loops (default=5) is the maximal number of loops carried out in the estimation procedure. (Each loops involves 
#   various methods to determine the high-likelihood region.)
# - until.convergence (default=T) is a Boolean indicating whether the estimation process should be terminated if 
#   there had been no improvement concerning the value of the target function between two consecutive loops. 
#   Otherwise, the estimation procedure is terminated according to the parameter "loops".
# - If print.output==T (default=T), interim results of the grid search and numerical optimization are printed
#   into the console throughout the estimation procedure.
# - If show.plots==T (default=T), profile log-likelihood plots are displayed at the end of the estimation procedure.
# - plot.title is the title of each plot if show.plots==T
# - pdf.file is an optional filename. If this is not missing and show.plots==T, the profile log-likelihoods will 
#   be plotted into this file.
# - If use.constraints=T, constraints on the densities of the populations will be applied.
# - "subgroups" is a list of sets of gene numbers. This parameter should be given only when the present call of 
#   stochprof.loop is based on a subanalysis of the subgroups of genes with non-fixed mu. The parameter is used 
#   only for calculation of the adjusted BIC which takes into account the number of parameters that had to be 
#   estimated during the whole estimation procedure: First, for each of the subclusters, and then for the final analysis.
     
     
   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)        
   get.range <- NULL
   stochprof.search <- NULL
   stochprof.results <- NULL
   penalty.constraint <- NULL
   calculate.ci <- NULL     
   rm(get.range)
   rm(stochprof.search)
   rm(stochprof.results)
   rm(penalty.constraint)
   rm(calculate.ci)     
     
   # check for obvious mistakes
   
   # model has to be defined
   if (!(model %in% c("LN-LN","rLN-LN","EXP-LN"))) {
      stop("stochprof.loop: unknown model.")
   }
   set.model.functions(model)   
   
   # dataset should be a matrix
   if (!is.matrix(dataset)) {
      stop("stochprof.loop: dataset has to be a  matrix.")
   }
   m <- ncol(dataset)
                  
   # check correct dimension of prev.result   
   if (!is.null(prev.result)) {
      if (model=="LN-LN") {
         correct.dim.pr <- (m+1)*TY+1
      }
      else if (model=="rLN-LN") {
         correct.dim.pr <- (m+2)*TY
      }
      else if (model=="EXP-LN") {
         if (TY==1) {
            correct.dim.pr <- m
         }
         else {
            correct.dim.pr <- (m+1)*TY+1
         }
      }   
   
      if (ncol(prev.result )!=correct.dim.pr) {
         stop("stochprof.loop: prev.result does not have the correct number of columns.")
      }
   }
   
   # check correct dimension of fixed.mu   
   if (!missing(fixed.mu)) {
      if (model %in% c("LN-LN","rLN-LN")) {
         correct.dim.fm <- TY*m
      }
      else if (model=="EXP-LN") {
         correct.dim.fm <- (TY-1)*m
      }    
   
      if (length(fixed.mu)!=correct.dim.fm) {
         stop("stochprof.loop: fixed.mu does not have the correct length.")
      }
   }   
     
   # if until.convergence==T, this Boolean indicates whether to exit the loop
   continue <- T
   # the currently best (i.e. minimal) value of the target function
   min.target <- 10^7
   # this variable will store the parameter at which the target function is evaluated and 
   # the value of the target function
   result <- prev.result
       
     
   for (i in 1:loops) {
      if (continue) {
         # try different methods to determine the high-likelihood region
         for (choice in c("none","mlci","quant","best")) {                
            # alternately carry out grid search and Nelder-Mead
            for (method in c("grid","optim")) {
               # determine appropriate parameter range ("high-likelihood region")
               this.range <- get.range(method=choice,prev.result=result,dataset=dataset,n=n,TY=TY,fix.mu=fix.mu,fixed.mu=fixed.mu)    
               if (is.null(this.range)) {
                  this.range <- par.range
               }         
               if (choice=="best") {
                  M <- 1
               }
               else if (method=="grid") {
                  M <- round(10^3/TY)
               }
               else {
                  M <- round(10/TY)
               } 

               # compute log-likelihood
               result <- stochprof.search(dataset=dataset,genenames=genenames,n=n,method=method,TY=TY,M=M,print.output=print.output,par.range=this.range,fix.mu=fix.mu,fixed.mu=fixed.mu,prev.result=result,use.constraints=use.constraints)
            }
         }
      }

      # has the algorithm converged?
      if (until.convergence) {
         result <- stochprof.results(prev.result=result,TY=TY,show.plots=F)
         if (!is.null(result)) {         
            new.min.target <- result[1,ncol(result)]
            if (round(new.min.target-min.target,4)==0) {
               continue <- F
            }
            else {
               min.target <- new.min.target
            }
         }
      }      
   }
         
   ################
   # final result #
   ################
   
   # log-likelihood for different parameter values
   pargrid <- stochprof.results(prev.result=result,TY=TY,show.plots=show.plots,plot.title=plot.title,pdf.file=pdf.file,fix.mu=fix.mu)
   if (is.null(pargrid)) {
      print("stochprof.loop: MLE not found.")
      return(invisible(NULL))
   }
   
   # maximum likelihood estimate
   mle <- pargrid[1,-ncol(pargrid)]
   # negative log-likelihood function
   loglikeli <- pargrid[1,ncol(pargrid)]
   names(loglikeli) <- NULL
   if (TY>1) {
      # penalty for constraint
      pen <- penalty.constraint(dataset,mle)  
   }
   else {
      pen <- 0
   }
   # BIC
   bic <- 2*loglikeli+log(length(dataset))*length(mle)
   # adjusted BIC
   if ((!missing(subgroups)) && (!is.null(subgroups))) {
      par.numb <- 0
      for (j in 1:length(subgroups)) {
         this.set <- subgroups[[j]]
         # number of parameters estimated for this subgroup
         if (model %in% c("LN-LN","EXP-LN")) {
            this.number <- TY*(length(this.set)+1)
         }
         else if (model=="rLN-LN") {
            this.number <- TY*length(this.set) + 3
         }
         par.numb <- par.numb + this.number
      }
      # plus number of parameters estimated in final analysis
      if (model=="LN-LN") {
         # F and sigma
         par.numb <- par.numb + 2
      }
      if (model=="rLN-LN") {
         # F and sigma=(sigma_1,...,sigma_TY)
         par.numb <- par.numb + 1 + TY
      }
      if (model=="EXP-LN") {
         # F and sigma and lambda=(lambda^(1),...,lambda^(m))
         par.numb <- par.numb + m + 2
      }
      adj.bic <- 2*loglikeli+log(length(dataset))*par.numb
   }
   else {
      adj.bic <- NULL
   }
   # confidence intervals
   ci <- calculate.ci(alpha=0.05,parameter=mle,dataset=dataset,n=n,TY=TY,fix.mu=fix.mu,fixed.mu=fixed.mu) 
   if (!is.null(ci)) {
      rownames(ci) <- colnames(pargrid)[-ncol(pargrid)]
   }
   
   cat("\n")
   cat("Maximum likelihood estimate (MLE):\n")
   print(mle)
   cat("\n")        
   cat("Value of negative log-likelihood function at MLE:\n")
   cat(loglikeli,"\n\n")
   cat("Violation of constraints:\n")
   if (pen==0) { cat("none\n\n") }
   else { cat(pen,"\n\n") }
   cat("BIC:\n")
   cat(bic,"\n\n")
   if ((!missing(subgroups)) && (!is.null(subgroups))) {
      cat("adjusted BIC:\n")
      cat(adj.bic,"\n\n")  
   }   
   cat("Approx. 95% confidence intervals for MLE:\n")
   print(ci)
   cat("\n")   
   cat("Top parameter combinations:\n")
   print(head(pargrid))
   cat("\n")
   
   final <- list(mle,loglikeli,ci,pargrid,bic,adj.bic,pen)
   names(final) <- c("mle","loglikeli","ci","pargrid","bic","adj.bic","pen")
   return(invisible(final))
}
