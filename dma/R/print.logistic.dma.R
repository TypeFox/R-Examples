print.logistic.dma <-
  function(x, ...){
    cat("Models:\n") 
    mattmp<-x$models
    ltmp<-"Model 1"
    if(nrow(mattmp)>1){
      for(i in 2:nrow(mattmp)){ltmp<-c(ltmp,paste("Model",i))}}
    rownames(mattmp)<-ltmp
    ltmp<-"Coef 1"
    if(ncol(mattmp)>1){
      for(i in 2:ncol(mattmp)){ltmp<-c(ltmp,paste("Coef",i))}}
    colnames(mattmp)<-ltmp
    print(mattmp) 
    cat("\nPosterior model probabilities:\n") 
    mattmp<-t(x$pmp)
    ltmp<-"Model 1"
    if(ncol(mattmp)>1){
      for(i in 2:ncol(mattmp)){ltmp<-c(ltmp,paste("Model",i))}}
    colnames(mattmp)<-ltmp
    print(mattmp)
  }
