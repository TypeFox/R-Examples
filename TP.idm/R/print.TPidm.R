print.TPidm <-
function(x, ...){
    if(!inherits(x,"TPidm")) stop("'x' must be of class 'TPidm'.")
    
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Parameters:","\n")
    cat("s=", x$s, "\n")
    cat("t=", x$t, "\n")
    cat("Method=", x$method, "\n")
    cat("CI=", x$CI, "\n")
    cat("CI transformation=", x$ci.transformation, "\n")
    cat(paste("Possible transitions:\n", sep=""))
    print(x$p.trans, row.names = FALSE)
    
    
    cat("\n")
    if(x$CI==TRUE){
      if(x$s==0){
        cat("Occupation probabilities at time t:", "\n", "\n")
        
        idm.info<-data.frame(transition=rownames(x$probs),
                             probs=x$probs[,1],
                             lower=x$probs[,2],
                             upper=x$probs[,3],
                             variance=x$probs[,4])
        
      }else{
        cat("Transition probabilities at time t:", "\n", "\n")
        
        idm.info<-data.frame(transition=rownames(x$probs),
                             probs=x$probs[,1],
                             lower=x$probs[,2],
                             upper=x$probs[,3],
                             variance=x$probs[,4])
        
      }
    }else{
      if(x$s==0){
        cat("Occupation probabilities at time t:", "\n", "\n")
        
        estimations<-c()
        for(j in 1:length(x$p.trans)){
          idx<-unlist(strsplit(x$p.trans[j], " "))
          estimations[j]<-x$all.probs[length(x$times),,j]
        }
        
        idm.info<-data.frame(transition=x$p.trans,
                             probs=estimations)
        
      }else{
        cat("Transition probabilities at time t:", "\n", "\n")
        
        estimations<-c()
        for(j in 1:length(x$p.trans)){
          idx<-unlist(strsplit(x$p.trans[j], " "))
          estimations[j]<-x$all.probs[length(x$times),,j]
        }
        
        idm.info<-data.frame(transition=x$p.trans,
                             probs=estimations)
        
      }
    }
    
    
    cat("\n")
    print(idm.info, row.names=FALSE)
    
    
  }
