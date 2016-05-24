boot.transitions <- function(transitions, iterations, by.stage.counts=FALSE, ...)
{
   ## check orderd fate, stage, and one or more fertility columns?
   t   <- iterations
   mat <- vector("list", t)          ## initialize a list to store matrices
   vec <- vector("list", t)          ## and stage vectors 
   lam <-numeric(t)                  ## and a vector for lambdas
   for (i in 1:t){ 
      if(by.stage.counts){
          ## create new data frame with resampled transitions by counts in original class vector
         boot <- do.call(rbind,                           
          lapply(split(transitions, transitions$stage, drop=TRUE),  
          function(x) x[sample(nrow(x), replace=TRUE),]))    
      } else{ 
         boot <- transitions[sample(nrow(transitions), replace=TRUE), ]
      }
      A <- projection.matrix(boot, ...) 
      # Nov 2011 - fixed bug noted by Tristan Lemke - if some stages are NA, they will be included in the summary below
      # and return an error when creating final list using rownames(A).  
      # vec[[i]] <- summary(boot$stage)  
      vec[[i]] <- table(boot$stage) 
      mat[[i]] <- as.vector(A)
      lam[i]   <- lambda(A)        
   }
   n <- dim(A)[1]
   boot.stage <- list(
     lambda= lam, 
     matrix= matrix(unlist(mat), byrow=TRUE, nrow=t,
          dimnames=list(1:t, paste("a", 1:n, rep(1:n,each=n), sep=""))),               
     vector= matrix(unlist(vec), byrow=TRUE, nrow=t, dimnames=list(1:t, names(vec[[1]])   ))  
   )
   boot.stage
}
