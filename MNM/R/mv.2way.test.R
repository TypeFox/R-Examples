###
###
###

`mv.2way.test` <- function(x, block, treatment,
                              score=c("identity","sign","rank"),
                              stand=c("outer","inner"),
                              method=c("approximation","permutation"),
                              n.simu=1000,
                              eps=1.0e-10,
                              n.iter=10000,
                              na.action=na.fail)
{
  DNAME=paste(deparse(substitute(x)),"by",deparse(substitute(treatment)),
    "within", deparse(substitute(block)))

  x<-na.action(x)
  if(!all(sapply(x, is.numeric))) stop("'x' must be numeric")
  x<-as.matrix(x)

  block<-na.action(block)
  if (!is.factor(block)) stop("'block' must be a factor")
  
  treatment<-na.action(treatment)
  if (!is.factor(treatment)) stop("'treatment' must be a factor")
  
  if(!all(table(block,treatment)==1))
    stop("Every treatment should be allocated once to every block")

  d <- dim(x)[2]
  if (d<2) stop("'x' must be at least bivariate")
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]

  #Sort the observations (permutation test would fail if
  #the data set was not sorted) 
  ordblock<-order(block)
  block<-block[ordblock]
  treatment<-treatment[ordblock]
  x<-x[ordblock,]
  
  score <- match.arg(score)
  stand <- match.arg(stand)
  method <- match.arg(method)

  res1<-switch(score,
         "identity"=
         {
           manova.identity(x=x,block=block,treatment=treatment,
                           method=method,nsim=n.simu)
         },
         "sign"=
         {
           manova.sign(x=x,block=block,treatment=treatment,
                       stand=stand,method=method,eps=eps,
                       maxiter=n.iter,nsim=n.simu)           
         },
         "rank"=
         {
           manova.rank(x=x,block=block,treatment=treatment,
                       stand=stand,method=method,eps=eps,
                       maxiter=n.iter,nsim=n.simu)
         }
         )

  NVAL<-paste("c(",paste(rep(0,d),collapse=","),")",sep="")
  names(NVAL)<-"location difference between some groups"
  ALTERNATIVE <- "two.sided"
    
  res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))

  class(res) <- "htest"    

  return(res)
}
