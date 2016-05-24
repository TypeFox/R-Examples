SurvtoLR <-
function(x){
    ## change Surv object to data.frame with L and R columns
    type<-attr(x,"type")
    if (type=="interval"){
        L<-R<-x[,1]
        R[x[,3]==0]<-NA
        L[x[,3]==2]<-0
        R[x[,3]==3]<-x[x[,3]==3,2]
    } else { stop(paste("Surv obj type='",type,"' unrecognized",sep=""))
    }
    out<-cbind(L,R)
    return(out)
  }
