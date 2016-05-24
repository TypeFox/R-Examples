analyze.toycluster <-
function(model="LN-LN",data.model="LN-LN",TY=2,preanalyze=T,show.plots=T,use.constraints=F) {
# Estimates the chosen model for the synthetic datasets provided in this package.
#
# - model is the model that is assumed for the data; can be "LN-LN", "rLN-LN" or "EXP-LN"
# - data.model is the model that was used for data simulation; can be "LN-LN", "rLN-LN" or "EXP-LN"
# - TY is the number of types of cells that is assumed in the model
# - if preanalyze=TRUE, the single-gene preanalysis as described below is carried out
# - if show.plots=TRUE, interim results are graphically displayed. This requires the user to confirm each new plot.

   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   toycluster.LNLN <- NULL
   toycluster.rLNLN <- NULL
   toycluster.EXPLN <- NULL
   get.range <- NULL
   rm(toycluster.LNLN)
   rm(toycluster.rLNLN)
   rm(toycluster.EXPLN)
   rm(get.range)
  
   # model and data model has to be defined
   if (!(model %in% c("LN-LN","rLN-LN","EXP-LN"))) {
      stop("Unknown model.")
   }
   if (!(data.model %in% c("LN-LN","rLN-LN","EXP-LN"))) {
      stop("Unknown model.")
   }    
   set.model.functions(model)
   
   ############
   # settings #
   ############
   
   genes <- 1:12   
   m <- length(genes)

   # length of vector mu   
   if (model=="LN-LN") {
      mu.length <- TY*m
   }
   else if (model=="rLN-LN") {
      mu.length <- TY*m
   }
   else if (model=="EXP-LN") {
      if (TY>1) {
         mu.length <- (TY-1)*m
      }
      else {
         mu.length <- 0
      }
   }   

   if (data.model=="LN-LN") {
      dataset <- toycluster.LNLN
   
      # Divide the 12-gene cluster into 3 groups of size 4, according to their order 
      # after the hierarchical clustering
      set.set <- list(c(7,5,2,8),c(1,3,4,10),c(9,6,12,11))      
   }
   else if (data.model=="rLN-LN") {
      dataset <- toycluster.rLNLN
   
      # Divide the 12-gene cluster into 3 groups of size 4, according to their order 
      # after the hierarchical clustering     
      set.set <- list(c(12,9,6,11),c(4,10,5,3),c(1,7,8,2)) 
   }
   else if (data.model=="EXP-LN") {
      dataset <- toycluster.EXPLN
   
      # Divide the 12-gene cluster into 3 groups of size 4, according to their order 
      # after the hierarchical clustering    
      set.set <- list(c(11,1,10,9),c(3,5,8,7),c(4,2,12,6))      
   }

   # according sets of indices with respect to vector "genes"   
   index.list <- list(length=length(set.set))
   for (j in 1:length(set.set)) {
      this.set <- set.set[[j]]
      k.vector <- NULL
      for (i in 1:length(this.set)) {
         this.gene <- this.set[i]   
         k <- which(genes==this.gene)
         k.vector <- c(k.vector,k)
      }
      index.list[[j]] <- k.vector
   }   
      
   #################
   # optional step #
   #################
   # Consider all genes separately first in order to get a rough idea about 
   # the location of the parameters. This might speed up the analysis of the larger clusters.       
   if ((preanalyze) && (mu.length>0)) {
      cat("###################################################\n")
      cat("## Preanalysis: Consider all genes individually. ##\n")
      cat("###################################################\n")      
      # number of parameters
      if (model=="LN-LN") {
         npar <- TY*2
         mu.indices <- TY:(npar-1)
         mu.names <- paste("mu_",1:TY,sep="")
      }
      else if (model=="rLN-LN") {
         npar <- TY*3-1
         mu.indices <- TY:(2*TY-1)
         mu.names <- paste("mu_",1:TY,sep="")
      }
      else if (model=="EXP-LN") {
         if (TY>1) {
            npar <- TY*2
            mu.indices <- TY:(2*TY-2)
            mu.names <- paste("mu_",1:(TY-1),sep="")
         }
         else {
            npar <- 1
            mu.indices <- NULL
            mu.names <- NULL
         }
      }
      
      # store the results for the single genes in this matrix:
      # lower and upper bounds of confidence intervals
      single.gene.lower <- matrix(nrow=npar,ncol=length(genes))
      single.gene.upper <- matrix(nrow=npar,ncol=length(genes))
                  
      # single gene estimation                  
      for (i in 1:length(genes)) {
         this.gene <- genes[i]
         
         cat("----------\n")
         cat(paste("| Gene",this.gene,"|\n"))
         cat("----------\n")

         single.res <- stochprof.loop(model=model,dataset=dataset[,i,drop=F],n=10,TY=TY,genenames=this.gene,fix.mu=F,loops=3,until.convergence=T,print.output=F,show.plots=show.plots,plot.title=paste("Gene",this.gene),use.constraints=use.constraints)
         
         if (!(is.null(single.res$ci))) {
            single.gene.lower[,i] <- (single.res$ci)[,1]
            single.gene.upper[,i] <- (single.res$ci)[,2]
         }
         else {
            other.interval <- get.range(method="quant",prev.result=single.res$pargrid,TY=TY,fix.mu=F)
            single.gene.lower[,i] <- other.interval[,1]
            single.gene.upper[,i] <- other.interval[,2]            
         }
      }
      
      colnames(single.gene.lower) <- paste("gene",genes)
      cn <- colnames(single.res$pargrid)
      rownames(single.gene.lower) <- cn[-length(cn)]
      if (length(mu.indices)>0) {
         rownames(single.gene.lower)[mu.indices] <- mu.names
      }
      if (model=="EXP-LN") {
         rownames(single.gene.lower)[nrow(single.gene.lower)] <- "lambda"
      }
      dimnames(single.gene.upper) <- dimnames(single.gene.lower)
            
      cat("--------------------------\n")
      cat("| Result of preanalysis: |\n")
      cat("--------------------------\n")      
      cat("lower bounds:\n")
      print(single.gene.lower)
      cat("\n")      
      cat("upper bounds:\n")
      print(single.gene.upper)
   }
      
   #############
   # main step #
   #############

   if (mu.length>0) {
      cat("#################################################\n")
      cat("## Main analysis: Consider subgroups of genes. ##\n")
      cat("#################################################\n")      
      
      # for each of the subgroups, estimate all parameters
      result.set <- list(length=length(set.set))
      
      for (set.index in 1:length(set.set)) {
         # consider the subsets of genes
         this.set <- set.set[[set.index]]
         # size of this subset
         this.m <- length(this.set)
         
         
         # strings for printing the results
         label <- "| Genes"
         for (v in 1:length(this.set)) {
            label <- paste(label,this.set[v])
         }
         label <- paste(label,"|")
         label2 <- NULL
         for (v in 1:nchar(label)) {
            label2 <- paste(label2,"-",sep="")
         }
         cat(label2,"\n")      
         cat(label,"\n")
         cat(label2,"\n")      
         
         if (preanalyze) {
            # range from which parameters should be drawn
            if (model=="LN-LN") {
               this.npar <- TY*(this.m+1)
               sigma.indices <- this.npar
               single.sigma.indices <- nrow(single.gene.lower)
               single.mu.indices <- TY:(nrow(single.gene.lower)-1)
            }
            else if (model=="rLN-LN") {
               this.npar <- TY*(this.m+2)-1
               sigma.indices <- ((this.m+1)*TY):((this.m+2)*TY-1) 
               single.sigma.indices <- (2*TY):(3*TY-1)  
               single.mu.indices <- TY:(2*TY-1)
            }
            else if (model=="EXP-LN") {
               if (TY>1) {
                  this.npar <- TY*(this.m+1)
                  sigma.indices <- (this.m+1)*(TY-1)+1        
                  single.sigma.indices <- 2*(TY-1)+1 
                  single.mu.indices <- TY:(2*(TY-1))
               }
               else {
                  this.npar <- this.m
                  sigma.indices <- NULL        
                  single.sigma.indices <- NULL  
                  single.mu.indices <- NULL          
               }
            }
            par.range <- matrix(nrow=this.npar,ncol=2)
            # initialize bounds for p
            if (TY>1) {
               par.range[1:(TY-1),1] <- 1
               par.range[1:(TY-1),2] <- 0
            }
            # initialize bounds for sigma  
            if (length(sigma.indices)>0) {    
               par.range[sigma.indices,1] <- 10^7
               par.range[sigma.indices,2] <- 0
            }
            if (model=="EXP-LN") {
               # initialize bounds for lambda          
               if (TY>1) {   
                  par.range[(this.m+1)*(TY-1)+1+(1:this.m),1] <- 10^7
                  par.range[(this.m+1)*(TY-1)+1+(1:this.m),2] <- 0         
               }
               else {
                  par.range[,1] <- 10^7
                  par.range[,2] <- 0                  
               }
            }
            # derive values from preanalysis
            for (i in 1:length(this.set)) {
               # consider each gene in the subset
               this.gene <- this.set[i]   
               # position of this gene in the overall dataset
               k <- which(genes==this.gene)
               # p: union of all confidence intervals for p in this subset
               if (TY>1) {
                  par.range[1:(TY-1),1] <- pmin(par.range[1:(TY-1),1],single.gene.lower[1:(TY-1),k])
                  par.range[1:(TY-1),2] <- pmax(par.range[1:(TY-1),2],single.gene.upper[1:(TY-1),k])            
               }
               # sigma: union of all confidence intervals for sigma in this subset
               if (length(sigma.indices)>0) {
                  par.range[sigma.indices,1] <- min(par.range[sigma.indices,1],single.gene.lower[single.sigma.indices,k])
                  par.range[sigma.indices,2] <- max(par.range[sigma.indices,2],single.gene.upper[single.sigma.indices,k])            
               }
               # mu: equal to confidence interval for mu for the respective gene           
               if (model %in% c("LN-LN","rLN-LN")) {
                  mu.indices <- TY-1+(0:(TY-1))*this.m+i
               }
               else if (model=="EXP-LN") {
                  if (TY>1) {
                     mu.indices <- TY-1+(0:(TY-2))*this.m+i 
                  }
                  else {
                     mu.indices <- NULL
                  }
               }                     
               if (length(mu.indices)>0) {
                  par.range[mu.indices,1] <- single.gene.lower[single.mu.indices,k]
                  par.range[mu.indices,2] <- single.gene.upper[single.mu.indices,k]
               }
               
               if (model=="EXP-LN") {
                  # lambda            
                  if (TY>1) {   
                     par.range[(this.m+1)*(TY-1)+1+i,1] <- pmin(par.range[(this.m+1)*(TY-1)+1+i,1],single.gene.lower[nrow(single.gene.lower),k])
                     par.range[(this.m+1)*(TY-1)+1+i,2] <- pmax(par.range[(this.m+1)*(TY-1)+1+i,2],single.gene.upper[nrow(single.gene.lower),k])         
                  }
                  else {
                     par.range[i,1] <- pmin(par.range[i,1],single.gene.lower[nrow(single.gene.lower),k])
                     par.range[i,2] <- pmax(par.range[i,2],single.gene.upper[nrow(single.gene.lower),k])       
                  }
               }
            }
         }
         else {
            par.range <- NULL
         }
         
         # set of indices
         indices <- index.list[[set.index]]
              
         # estimation for entire subgroup
         this.result <- NULL      
         while (is.null(this.result)) {
            this.result <- stochprof.loop(model=model,dataset=dataset[,indices,drop=F],n=10,TY=TY,par.range=par.range,genenames=this.set,fix.mu=F,loops=10,until.convergence=T,print.output=F,show.plots=show.plots,plot.title=paste("Subgroups of genes: group number",set.index),use.constraints=use.constraints)
         }
         result.set[[set.index]] <- this.result
      }
   }
   
   
   ##############
   # final step #
   ##############
   
   cat("################################################################\n")
   cat("## Final analysis: Fix mu, estimate p and sigma (and lambda). ##\n")
   cat("################################################################\n")         

   # Fix mu to the values determined in the subclusters.
   # There remain the following parameters to estimate: 
   # p and sigma in the LN-LN and rLN-LN model, and
   # p, sigma and lambda in the EXP-LN model.
   if (mu.length>0) { 
      fixed.mu <- vector(length=mu.length)
      for (i in 1:length(result.set)) {
         # result of this subgroup
         this.result <- result.set[[i]]
         this.mle <- this.result$mle
         # indices of genes in this subgroup
         indices <- index.list[[i]]
         # according indices in fixed.mu of entire cluster
         if (model=="LN-LN") {
            this.mu.indices <- TY:(length(this.mle)-1)
         }
         else if (model=="rLN-LN") {
            this.mu.indices <- TY:(length(this.mle)-TY)
         }
         else if (model=="EXP-LN") {
            if (TY>1) {
               this.mu.indices <- TY:(length(this.mle)-1-this.m)
            }
            else {
               this.mu.indices <- NULL         
            }
         }
         fm.indices <- NULL
         for (j in 1:(mu.length/m)) {
            fm.indices <- c(fm.indices,(j-1)*m+indices)
         }
         fixed.mu[fm.indices] <- this.mle[this.mu.indices]
      }      
      # estimation of p, sigma and lambda for mu fixed to the values estimated in the subgroups
      result <- stochprof.loop(model=model,dataset=dataset[,genes,drop=F],n=10,TY=TY,genenames=genes,fix.mu=T,fixed.mu=fixed.mu,loops=10,until.convergence=T,print.output=F,show.plots=show.plots,plot.title="Entire Cluster",use.constraints=use.constraints,subgroups=set.set)
   }   
   else {
      # estimation of p and lambda (there is no mu)
      result <- stochprof.loop(model=model,dataset=dataset[,genes,drop=F],n=10,TY=TY,genenames=genes,fix.mu=F,loops=10,until.convergence=F,print.output=F,show.plots=F,plot.title="Entire Cluster",use.constraints=use.constraints)
   }      
   
   
   return(invisible(result))
}
