pop.projection <- function(A,n,iterations=20)
{
   x<-length(n)
   t <- iterations  
   stage <- matrix(numeric(x*t), nrow=x)          ## initialize a matrix to store projected stage vectors
   pop <- numeric(t)                              ## and numeric vectors for projected population size
   change <- numeric(t-1)                         ## and proportional changes in pop sizes

   for( i in 1:t)
   {
      stage[,i] <- n  
      pop[i]   <- sum(n)          
      if(i>1){change[i-1] <- pop[i]/pop[i-1] }    ## calculates proportional changes in pop size
      n <- A %*% n                                ## multiply matrix A by size vector n and 
   }                                              ## set n equal to new vector
 

   rownames(stage)<-rownames(A)                   ## and add row names.
   colnames(stage)<-0:(t-1)                       ## start counting at zero?
    w <- stage[,t]                                ## Estimate stable stage from last size iteration

   pop.proj <- list(  
     lambda = pop[t]/pop[t-1],                    ## Change between the LAST two population counts  
     stable.stage  = w/sum(w),
     stage.vectors = stage,   
     pop.sizes = pop,
     pop.changes = change 
   )
   pop.proj
}

