stochprof.search.EXPLN <-
function(dataset,n,TY,method="grid",M=10,par.range=NULL,prev.result=NULL,fix.mu=F,fixed.mu,genenames=NULL,print.output=F,use.constraints=F) {
# Calculates the log likelihood function of all model parameters for a given dataset
# at certain parameter values. The so-obtained values are returned in a matrix with
# the following entries: Each row corresponds to one parameter combination. All columns 
# but the last one contain the parameter values at which the log likelihood function has 
# been computed. The column names are the parameter names. The last column ("target") is the 
# negative log likelihood function computed at the respective parameter vector. For numerical 
# reasons, this target value is set to the minimum of 10^7 and the actual value.
# 
# The values at which the target function is calculated are randomly drawn from some range 
# specified by "par.range". If method=="grid", the target function is simply evaluated
# at such a randomly drawn parameter vector. If method=="optim", this randomly drawn vector is
# passed to the Nelder-Mead algorithm as a starting value in order to search for a local 
# maximum around it.
#
# Parameters:
#
# - dataset is a matrix which contains the cumulated expression data over all cells in a tissue sample.
#   Columns represent different genes, rows represent different tissue samples.
# - n is the number of cells taken from each tissue sample.
# - TY is the number of types of cells that is assumed in the stochastic model.
# - method (default="grid") determines whether a grid search or the Nelder-Mead algorithm should be applied: 
#   If method=="grid", the log likelihood function is simply evaluated at certain parameter values that are
#   randomly drawn.
#   If method=="optim", a Nelder-Mead search starts at a randomly drawn set of parameter values in order to 
#   find a local maximum. The resulting locally optimal parameter is stored in the results matrix as one row.
# - M (default=10) is the number of randomly drawn parameter combinations.
# - par.range (default=NULL) is the range from which the parameter values should be randomly drawn. This is 
#   a matrix with the number of rows being equal to the number of model parameters. The first column contains 
#   the lower bound, the second column the upper bound. If par.range==NULL, some rather large range is defined.
# - prev.result (default=NULL) can contain results from former calls of this function.
# - fix.mu (default=F) indicates whether the log-means are kept fixed in the estimation procedure or whether
#   they are to be estimated.
# - fixed.mu (no default, needs to be specified only when fix.mu==T) is a vector containing the values to which
#   the log-means should be fixed if fix.mu==T. The order of components is as follows:
#   (mu_type_1_gene_1, mu_type_1_gene_2, ..., mu_type_2_gene_1, mu_type__gene_2, ...)
# - genenames (default=NULL) are the names of the genes in the dataset.
#   For genenames==NULL, the genes will simply be enumerated according to the column numbers in the dataset.
# - If print.output==T (default=F), interim results of the grid search and numerical optimization are printed
#   into the console throughout the estimation procedure.
# - If use.constraints==T, constraints on the densities of the populations will be applied.

   
   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)        
   d.sum.of.mixtures <- NULL
   backtransform.par <- NULL
   penalty.constraint <- NULL
   draw.parameters <- NULL
   transform.par <- NULL
   rm(d.sum.of.mixtures)
   rm(backtransform.par)
   rm(penalty.constraint)
   rm(draw.parameters)
   rm(transform.par)
   
   
   if (M==0) return(NULL)  
      
   ####################
   # general settings #              
   ####################

   # number of genes
   m <- ncol(dataset)
   # gene names
   if (is.null(genenames)) {
      genenames <- 1:m
   }
   
   # names of variables
   if (TY==1) {
      varnames <- c(paste("lambda_",genenames,sep=""))
   }
   else {
      varnames <- paste("p_",1:(TY-1),sep="")
      for (i in 1:(TY-1)) {
         varnames <- c(varnames,paste("mu",i,"gene",genenames,sep="_"))
      }
      varnames <- c(varnames,"sigma")
      varnames <- c(varnames,paste("lambda_",genenames,sep=""))
   }
        
   ###########################################
   ## ML estimation: define target function ##
   ###########################################
   
   # loglikelihood for one gene
   loglikeli <- function(y,p,mu,sigma,lambda) {
      # p is of length TY
      # mu and sigma are of length TY-1 (so possibly NULL)
      # lambda is of length 1
      max(-10^7,sum(d.sum.of.mixtures(y,n,p,mu,sigma,lambda,logdens=T)))
   }
      
   # this function will be minimized (the function "to.minimize" below just 
   # changes the parameterisation)
   target.function <- function(p,mu,sigma,lambda) {  
      # p is of length TY
      # mu is of length (TY-1)*m
      # sigma is of length 1 (if TY>1) or NULL
      # lambda is of length m
      
      # build mu.matrix such that the g.th row contains the values for gene g
      if (!is.null(mu)) {
         mu.matrix <- matrix(mu,byrow=F,ncol=TY-1)
      }
      
      # consider negative log likelihood because the target function will be minimized
      this.sum <- 0
      for (g in 1:m) {
         if (TY>1) {
            this.sum <- this.sum - loglikeli(dataset[,g],p,mu.matrix[g,,drop=T],rep(sigma,TY-1),lambda[m])
         }
         else {
            this.sum <- this.sum - loglikeli(dataset[,g],p,NULL,NULL,lambda[m])
         }
      }
      return(this.sum)
   }
   
   # For identifiability purposes, we require the mu-values of the first gene to be
   # in descending order. Otherwise, the parameter combinations (p,mu1_gene1,mu2_gene1,sigma)
   # and (1-p,mu2_gene1,mu2_gene1,sigma) would yield identical values of the log-likelihood 
   # function.
   # The randomly drawn parameters will always fulfil the above requirement. However, during the
   # Nelder-Mead optimization procedure it might happen that it is violated. In that case, the 
   # optim algorithm might jump back and forth between two equivalent states.  
   # In order to avoid this, a penalty term is introduced which is added to the actual target 
   # function. This penalty is positive if mu_1^1 >= mu_2^1 >= ... >= mu_TY^1 is not fulfilled, 
   # and 0 otherwise.   
   penalty.mu <- function(mu,m,lambda=100) {
      # this function should only be called if TY>2
      
      # build a matrix such that the g.th column contains the mu values for gene g
      mu <- matrix(mu,byrow=T,ncol=m)
      # mu for gene 1
      mu.g1 <- mu[,1]    
      # penalty
      differences <- mu.g1[-1]-mu.g1[-length(mu.g1)]
      pen <- pmax(0,differences)
      return(lambda*sum(pen^2))
   }   
   
   # this function should be minimized
   to.minimize <- function(theta) {   
      # theta=(w_1,...,w_{T-1},mu,log(sigma),log(lambda)).
      # w_i is one-dimensional;
      # mu is multi-dimensional (for TY>1);
      # sigma is one-dimensional (for TY>1);
      # lambda is m-dimensional.
      # If fix.mu==F and TY>1, then mu is (TY-1)*m-dim., otherwise zero-dimensional.
            
      # backtransformation           
      # afterwards, back.theta is full-dim, incl. mu
      back.theta <- backtransform.par(this.par=theta,m=m,fix.mu=fix.mu,fixed.mu=fixed.mu)
      
      if (TY>1) {
         p <- back.theta[1:(TY-1)]
         p <- c(p,1-sum(p))

         mu <- back.theta[TY:((m+1)*(TY-1))]
         sigma <- back.theta[(m+1)*(TY-1)+1]

         lambda <- back.theta[(m+1)*(TY-1)+1+(1:m)]         
      }
      else {
         p <- 1
         mu <- NULL
         sigma <- NULL
         lambda <- back.theta
      }      
      
            
      # penalty
      pen.mu <- 0
      if (TY>2) {
         pen.mu <- penalty.mu(mu,m)
      }
      pen.constr <- 0
      if ((TY>1) && (use.constraints)) {
         pen.constr <- penalty.constraint(dataset,parameter=c(p[-TY],mu,sigma,lambda))
      }
            
      # target
      a <- target.function(p,mu,sigma,lambda) + pen.mu + pen.constr
      a <- min(10^7,a)      
      return(a)
   }
   
   ####################
   # previous results #
   ####################
   
   # are there previous results already?
   if (is.null(prev.result)) {
      # no, there aren't
      all.results <- matrix(nrow=0,ncol=length(varnames)+1)
      colnames(all.results) <- c(varnames,"target")
   }    
   else {
      # yes, there are
      all.results <- prev.result
   }
   
                
   ###################################
   # parameter ranges to be searched #
   ###################################

   if (is.null(par.range)) {
      # determine some range
      ranges <- matrix(NA,ncol=2,nrow=length(varnames))

      if (TY>1) {
         # p      
         ranges[1:(TY-1),1] <- 0
         ranges[1:(TY-1),2] <- 1
         # mu
         ranges[TY:((m+1)*(TY-1)),1] <- -4.5
         ranges[TY:((m+1)*(TY-1)),2] <- 2.5
         # sigma
         ranges[TY:((m+1)*(TY-1))+1,1] <- 0.01
         ranges[TY:((m+1)*(TY-1))+1,2] <- 1                    
      }
      # lambda
      ranges[nrow(ranges)-m+(1:m),1] <- 0.01
      ranges[nrow(ranges)-m+(1:m),2] <- 50     
   }
   else {
      ranges <- par.range
   }  
   if ((fix.mu) && (TY>1)) {
      ranges[TY:((m+1)*(TY-1)),1] <- fixed.mu
      ranges[TY:((m+1)*(TY-1)),2] <- fixed.mu
   }      
      
   ###################################
   ## estimation: optimization step ##
   ###################################        
       
   #--------------#
   # optim method #
   #--------------#   
   if (method=="optim") {
      for (i in 1:M) {
         # draw starting value
         par0 <- draw.parameters(ranges,m) # full-dim. 
         theta0 <- transform.par(this.par=par0,m=m,fix.mu=fix.mu) # lower-dim. if fix.mu==T

         theta0[theta0==-Inf] <- -10^7
         theta0[theta0==Inf] <- 10^7         

         # numerically optimize
         if (length(theta0)==1) {
            result <- optim(theta0,fn=to.minimize,control=list(maxit=10^5),method="Brent",hessian=F, lower=-10^7, upper=10^7)
         }
         else {
            result <- optim(theta0,fn=to.minimize,control=list(maxit=10^5),method="Nelder-Mead",hessian=F)         
         }
         
         # result
         this.theta <- result$par
         this.par <- backtransform.par(this.par=this.theta,m=m,fix.mu=fix.mu,fixed.mu=fixed.mu) # full-dim              
         this.value <- result$value        
         
         # attach new result to all former ones 
         all.results <- rbind(all.results,c(this.par,this.value))
         
         if (print.output) {
            cat("---\n")
            cat("Start optim at:\n")
            cat(par0,"\n")
            cat("Arrived at:\n")
            cat(this.par,"\n")
         }         
      }
   }
   #-------------#
   # grid search #
   #-------------#   
   else if (method=="grid") {                
      for (i in 1:M) {       
         # randomly draw parameter 
         this.par <- draw.parameters(ranges,m) # full-dim
         this.theta <- transform.par(this.par=this.par,m=m,fix.mu=fix.mu) # lower-dim if fix.mu==T

         if (print.output) {
            cat("---\n")
            cat("Compute grid at:\n")
            cat(this.par,"\n")      
         }
         
         # target function
         this.value <- to.minimize(this.theta)  
         # attach new result to all former ones        
         all.results <- rbind(all.results,c(this.par,this.value))
      }
   }
   
   return(all.results)
}
