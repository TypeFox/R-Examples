alphascore <- function(xobs, yobs, xeval=xobs, yeval=yobs, alpha=0.95)
{
 n1<-nrow(xobs)    # number of observations
 p.input<-ncol(xobs)     # number of inputs
 q.output<-ncol(yobs)   # number of outputs
 n2<-nrow(xeval)     # number of evaluation points

# verification of some conditions
if(n1!=nrow(yobs)) stop("xobs and yobs have not the same number of observations")
if(n2!=nrow(yeval)) stop("xeval and yeval have not the same number of observations")

if(length(alpha)==1) alpha<-rep(alpha,n2)
if(any(alpha<=0) | any(alpha>1)) stop("alpha must be included in the interval ]0;1]")

 # Call to the C.code
 res<-.C("orderalpha",as.integer(n1), as.integer(n2), as.integer(p.input), as.integer(q.output),
 as.double(t(xobs)), as.double(t(yobs)), as.double(t(xeval)), as.double(t(yeval)),
 res1=as.double(matrix(0,n2,1)),res2=as.double(matrix(0,n2,1)),
 res3=as.double(matrix(0,n2,1)), res4=as.double(matrix(0,n2,1)),
 res5=as.double(matrix(1,n2,1)),res6=as.double(matrix(0,n2,1)),
 res7=as.double(matrix(0,n1,1)),res8=as.double(matrix(1,n1,1)),
  res9=as.double(matrix(0,n1,1)), as.double(alpha) ,PACKAGE="frontiles")

 # Extraction of the results
 orderalpha.score<-as.data.frame(cbind(1/res$res1,res$res3,res$res5))
 colnames(orderalpha.score)<-c("output","input","hyper")
 
 orderalpha.score[orderalpha.score==-1]<-NA

return(orderalpha.score)
}

