superpc.train<-
  function (data, type=c("survival","regression"), s0.perc=NULL){
    
# computes feature scores for supervised pc analysis
    
  this.call <- match.call()
 type <- match.arg(type)

  
  if (is.null(data$censoring.status) & type=="survival") {
 stop("Error: survival specified but  censoring.status is null")
  }
  
 if (!is.null(data$censoring.status) & type=="regression") {
 stop("Error: regression specified but  censoring.status is  non-null")
  }

  if (type=="survival") {
junk<- coxfunc(data$x, data$y, data$censoring.status, s0.perc=s0.perc)
feature.scores<-junk$tt
   }
  else {
junk<- cor.func(data$x, data$y, s0.perc=s0.perc)
feature.scores<-junk$tt
  }

  
  junk <- list( feature.scores=feature.scores, 
               type=type, s0.perc=s0.perc, 
               call = this.call)

  
  class(junk) <- "superpc"
  return(junk)

}

