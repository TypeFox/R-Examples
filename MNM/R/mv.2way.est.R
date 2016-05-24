###
###
###

`mv.2way.est` <- function(x, block, treatment,
                          score=c("identity","sign","rank"),
                          stand=c("outer","inner"),
                          eps=1.0e-10,
                          n.iter=1000,
                          na.action=na.fail)
{
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

  #Sort the observations
  ordblock<-order(block,treatment)
  block<-block[ordblock]
  treatment<-treatment[ordblock]
  x<-x[ordblock,]
  
  score <- match.arg(score)
  stand <- match.arg(stand)
 
  res<-switch(score,
         "identity"=
         {
           estimate.identity(x=x,block=block,treatment=treatment)
         },
         "sign"=
         {
           estimate.sign(x=x,block=block,treatment=treatment,stand=stand,eps=eps,maxiter=n.iter)           
         },
         "rank"=
         {
           estimate.rank(x=x,block=block,treatment=treatment,stand=stand,eps=eps,maxiter=n.iter)
         }
         )

  class(res) <- "mvcloc"
  ll<-length(res)
  for(i in 1:ll){
    names(res)[i] <- paste("Component ",i," (",res[[i]]$dname,")",sep="")
  }
  return(res)
}
