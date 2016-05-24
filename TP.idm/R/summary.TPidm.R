summary.TPidm <-
function(object,...){
    if(!inherits(object,"TPidm")) stop("'object' must be of class 'TPidm'.")
    
    cat("Call:\n")
    print(object$call)
    cat("\n")
    cat("Parameters:","\n")
    cat("\n")
    cat(paste("Time s:\n", sep=""))
    print(object$s, row.names = FALSE)
    cat("\n")
    cat(paste("Time t:\n", sep=""))
    print(object$t, row.names = FALSE)
    cat("\n")
    cat(paste("Method:\n", sep=""))
    print(object$method, row.names = FALSE)
    cat("\n")
    cat(paste("CI:\n", sep=""))
    print(object$CI, row.names = FALSE)
    cat("\n")
    cat(paste("CI transformation:\n", sep=""))
    print(object$ci.transformation, row.names = FALSE)
    cat("\n")
    cat(paste("Possible transitions:\n", sep=""))
    print(object$p.trans, row.names = FALSE)
    
    cat("\n")
    if(object$CI==TRUE){
      if(object$s==0){
        cat("Occupation probabilities at time t:", "\n", "\n")
        
        idm.info<-data.frame(transition=rownames(object$probs),
                             probs=object$probs[,1],
                             lower=object$probs[,2],
                             upper=object$probs[,3],
                             variance=object$probs[,4])
        
      }else{
        cat("Transition probabilities from s to t:", "\n", "\n")
        
        idm.info<-data.frame(transition=rownames(object$probs),
                             probs=object$probs[,1],
                             lower=object$probs[,2],
                             upper=object$probs[,3],
                             variance=object$probs[,4])
        
      }
    }else{
      if(object$s==0){
        cat("Occupation probabilities at time t:", "\n", "\n")
        
        estimations<-c()
        for(j in 1:length(object$p.trans)){
          idx<-unlist(strsplit(object$p.trans[j], " "))
          estimations[j]<-object$all.probs[length(object$times),,j]
        }
        
        idm.info<-data.frame(transition=object$p.trans,
                             probs=estimations)
        
      }else{
        cat("Transition probabilities from s to t:", "\n", "\n")
        
        estimations<-c()
        for(j in 1:length(object$p.trans)){
          idx<-unlist(strsplit(object$p.trans[j], " "))
          estimations[j]<-object$all.probs[length(object$times),,j]
        }
        
        idm.info<-data.frame(transition=object$p.trans,
                             probs=estimations)
        
      }
    }
    
    cat("\n")
    print(idm.info, row.names=FALSE)
    
    
  }
