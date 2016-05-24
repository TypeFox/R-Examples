ordermscore <- function(xobs, yobs, xeval=xobs, yeval=yobs, m=30)
{
# initialization
 n1<-nrow(xobs)    # number of observations
 p.input<-ncol(xobs)     # number of inputs
 q.output<-ncol(yobs)   # number of outputs
 n2<-nrow(xeval)     # number of evaluation points
 

 # verification of some conditions
if(n1!=nrow(yobs)) stop("xobs and yobs have not the same number of observations")
if(n2!=nrow(yeval)) stop("xeval and yeval have not the same number of observations")

if(m<=0) stop("m must be positive")


 # Call to the C.code
 res<-.C("orderm",as.integer(n1), as.integer(n2), as.integer(p.input), as.integer(q.output),
 as.double(t(xobs)), as.double(t(yobs)), as.double(t(xeval)), as.double(t(yeval)),
 res1=as.double(matrix(0,n2,1)),
 res3=as.double(matrix(0,n2,1)),
 res5=as.double(matrix(0,n2,1)),
 res7=as.double(matrix(0,n1,1)),res8=as.double(matrix(1,n1,1)),
  res9=as.double(matrix(0,n1,1)), as.double(m) ,PACKAGE="frontiles")

   # Extraction of the results
 orderm.score<-as.data.frame(cbind(1/res$res1,res$res3,1/res$res5))
 colnames(orderm.score)<-c("output","input","hyper")

 orderm.score[orderm.score==-1]<-NA

  return(orderm.score)
}

