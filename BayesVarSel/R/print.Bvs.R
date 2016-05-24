print.Bvs <-
function(x, ...){
  cat("\n")
  cat("Call:\n")
  print(x$call)
  #cat("\nThis is the result for a model selection problem with ")
  #cat(x$p-1)
  pp<-2^x$p-1
  n.keep<-dim(x$modelsprob)[1]
  #cat(" covariates and ")
  #cat(x$n)
  #cat(" observations\n")
  #cat("The potential covariates are:\n")
  #cat(x$variables[-1])
  #if(!is.null(x$time)){
  #  cat("\nComputational time: ")
  #  cat(x$time)
  #  cat(" seconds.\n")
  #}
  if(x$method=="gibbs"){
    cat("\nAmong the visited models, the model with the largest probability contains: \n")
    print(names(which(x$HPMbin==1)))
  }
  if(x$method!="gibbs"){
    if(n.keep<=10){
      cat(paste("\nThe",n.keep,"most probable models and their probabilities are:\n",sep=" "))
      print(x$modelsprob)
    }
    if(n.keep>10){
      cat("\nThe 10 most probable models and their probabilities are:\n")
      print(x$modelsprob[1:10,])
      cat("\n(The remanining", n.keep-10, "models are kept but omitted in this print)")
    }
  }
   cat("\n")
    
  }
    
