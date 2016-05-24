`CInp` <-
function(x, ...){UseMethod("CInp")}

`CInp.default` <-
function(x, conf.level=0.95, alternative="two.sided", ...)
{
alternative <- match.arg(alternative, choices=c("two.sided","less","greater"))

args<-list(...)

DataMatrix <- x
N <- nrow(DataMatrix)
k <- round(conf.level*N,0)

switch(alternative,

"two.sided"={
probs<-c((1-conf.level)/2, 1-(1-conf.level)/2)
CIs <- t( apply( X=DataMatrix, MARGIN=2, 
 FUN=function(x){quantile(x=x, probs=probs)} ))
},

"less"={
probs<-c(conf.level)
upper <- t( apply( X=DataMatrix, MARGIN=2, 
 FUN=function(x){quantile(x=x, probs=probs)} ))
CIs<-cbind(-Inf, upper)
},

"greater"={
probs<-c(1-conf.level)
lower <- t( apply( X=DataMatrix, MARGIN=2, 
 FUN=function(x){quantile(x=x, probs=probs)} ))
CIs<-cbind(lower,Inf)
}
)
# end of switch

estimate <- apply(X=DataMatrix, MARGIN=2, median)

colnames(CIs)<-c("lower","upper")

out<-list(
conf.int=CIs,
estimate=estimate,
x=x,
k=k,
N=N,
conf.level=conf.level,
alternative=alternative)

class(out)<-"CInp"

return(out)

}


`CInp.CCRatio` <-
function(x,...)
{
args<-list(...)

args$x<-x$chains

out<-do.call("CInp.default", args)

out$x<-x

return(out)

}

`CInp.CCDiff` <-
function(x,...)
{
args<-list(...)

args$x<-x$chains

out<-do.call("CInp.default", args)

out$x<-x

return(out)

}

`CInp.bugs` <-
function(x, conf.level=0.95, alternative="two.sided", whichp=NULL, ...)
{

args<-list(...)

sl<-x$sims.list

if(is.null(whichp))
{
mat<-x$sims.matrix
}
else{
 namsl<-names(sl)
 if(!whichp %in% namsl)
  {stop("whichp could not be found in the parameter list of the openbugs object")}

  if(length(whichp)==1)
   {
   mat<-sl[[whichp]]
   }
  if(length(whichp)>1)
   {
   mat<-matrix(nrow=x$n.sims)
   for (i in seq(along.with=whichp))
    {
     mat<-cbind(mat,x$sims.list[[whichp[i]]])
    }
   }
 }

args$x<-mat
args$conf.level<-conf.level
args$alternative<-alternative

out<-do.call("CInp.default", args)

return(out)
}

