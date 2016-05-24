
test.Rhat <- function(samples,cutoff,params.omit,verbose=TRUE){
  
  params <- colnames(samples[[1]])
  expand <- sapply(strsplit(params, "\\["), "[", 1)
  
  gd <- function(hold){
    r <- try(gelman.diag(hold, autoburnin=FALSE)$psrf[1], silent=TRUE)
    if(inherits(r, "try-error") || !is.finite(r)) {
      r <- NA
    }
    return(r)
  }
  
  failure <- FALSE
  index <- 1
  while (failure==FALSE && index <= length(params)){
    
    if(!expand[index]%in%params.omit){
      test <- gd(samples[,index])
    } else {test <- 1}
    
    if(is.na(test)){test <- 1}
   
    if(test>cutoff){failure=TRUE
    } else {index <- index + 1}
  }
  
  if(failure==TRUE&verbose){
    cat('.......Convergence check failed for parameter \'',params[index],'\'\n',sep="")
  }
  if(failure==FALSE&verbose){
    cat('.......All parameters converged.','\n\n')
  }
  
  return(failure)
  
}