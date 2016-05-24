#For the subset of the simplex defined by
# A1 p = b1, A2 p <= b2 and A3 p >= b3
# where the Ai's are matrices and the bi's
#vectors of nonnegative real numbers this
#function uses the Metroplis-Hastings algorithm

constrppmn<-function(A1,A2,A3,b1,b2,b3,initsol,reps,ysamp,burnin)
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
     out<-polyaest(A1,A4,b4,initsol,reps,ysamp,burnin)
     return(out)
  }

