print.CEoptim<-function(x,...){
  
  if(class(x)!="CEoptim")
    stop("The object x is not the result from CEoptim function")
  
  args <- list(...)
  OutPut<- list(optimizer=TRUE, optimum=TRUE, termination=TRUE, states=FALSE,states.probs=FALSE)  
  
  if((!is.null(args$optimizer)&&args$optimizer==TRUE)||
       (!is.null(args$optimum)&&args$optimum==TRUE)||
       (!is.null(args$termination)&&args$termination==TRUE)){ 
    OutPut<- list(optimizer=FALSE, optimum=FALSE, termination=FALSE, states=FALSE,states.probs=FALSE)
  }  
  nameOutPut <- names(OutPut)   #name inner list
  nameInput <- names(args)
  if (length(noNms <- nameInput[!nameInput %in% nameOutPut])!=0) {
    stop("unknown parameter(s) output: ",
         paste(noNms, collapse = ", "))
  }  else OutPut[nameInput]<-args
  
  
  optimizer <- OutPut$optimizer
  optimum <- OutPut$optimum
  termination <- OutPut$termination
  states <- OutPut$states
  states.probs <- OutPut$states.probs
  
  
  if(optimizer==TRUE){ 
    if(length(x$optimizer$continuous)>0) cat("Optimizer for continuous part:","\n",x$optimizer$continuous,"\n")
    if(length(x$optimizer$discrete)>0) cat("Optimizer for discrete part:","\n",x$optimizer$discrete,"\n")
  } 
  
  if(optimum==TRUE){   
    cat("Optimum:","\n",x$optimum,"\n")
  }
  
  if(termination==TRUE){
    cat("Number of iterations:","\n",x$termination$niter,"\n")
    cat("Total number of function evaluations:","\n",x$termination$nfe,"\n")
    cat("Convergence:","\n",x$termination$convergence,"\n")
  }  
  
  
  if(states==TRUE){
    cat("states:","\n")
    print(x$states)
  }
  
  if(states.probs==TRUE){
    cat("Categorical sampling probabilities:","\n")
    print(x$states.probs)
    
  }
  
}
