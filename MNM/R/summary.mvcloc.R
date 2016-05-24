`summary.mvcloc` <-
function(object,..., digits=4)
    {     
    p<-length(object)
    d<-length(c(object[[1]]$location))
    A<-matrix(0,p,d)
    rnames<-NULL
    for(i in 1:p){
      A[i,]<-object[[i]]$location
      rnames[i]<-object[[i]]$dname
    }
    vcov<-object[[1]]$vcov
    rownames(A)<-rnames
    colnames(A)<-1:d
    cat("The", object[[1]]$est.name, "\n")
    print(round(A,digits))
    if(is.numeric(vcov)){
      cat("\n")
      cat("The covariance matrix")
      cat("\n")
      print(round(vcov,digits))
    }
  
    invisible(list(A,vcov))
    }
