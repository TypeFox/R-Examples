#Generates k probability vectors in "steps" of
#size step.

constrppprob<-function(A1,A2,A3,b1,b2,b3,initsol,step,k)
  {
    checkconstr(A1,A2,A3,b1,b2,b3)
    if(!is.null(A3)) {
      A4<-rbind(A2,-A3)
      b4<-c(b2,-b3)
    }
     else {
       A4<-A2
       b4<-b2
     }
    nrows<-nrow(A4)
    out<-probvect(A1,A4,nrow(A4),b4,initsol,step,k)
    return(out)
  }
       
