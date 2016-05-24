`SCSnp` <-
function(x, ...){UseMethod("SCSnp")}

`SCSnp.default` <-
function(x, conf.level=0.95, alternative="two.sided", ...)
{
alternative <- match.arg(alternative, choices=c("two.sided","less","greater"))

DataMatrix <- x
N <- nrow(DataMatrix)
k <- round(conf.level*N,0)
RankDat <- apply(DataMatrix,2,rank)

switch(alternative,

"two.sided"={
W1 <- apply(RankDat,1,max)
W2 <- N + 1 - apply(RankDat,1,min)

Wmat <- cbind(W1,W2)
w <- apply(Wmat,1,max)
tstar <- round(sort(w)[k],0)

SCI <- function(x)
{
 sortx <- sort(x)
 cbind(sortx[N+1-tstar],sortx[tstar])
}

SCS <- t(apply(DataMatrix,2,SCI))
},

"less"={
W1 <- apply(RankDat,1,max)
tstar <- round(sort(W1)[k],0)

SCI <- function(x)
{
 sortx <- sort(x)
 cbind(-Inf, sortx[tstar])
}

SCS<-t(apply(DataMatrix,2,SCI))
},

"greater"={
W2 <- N + 1 - apply(RankDat,1,min)
tstar <- round(sort(W2)[k],0)

SCI <- function(x)
{
 sortx <- sort(x)
 cbind(sortx[N+1-tstar], Inf)
}

SCS<-t(apply(DataMatrix,2,SCI))

}
)
# end of switch

estimate<-apply(DataMatrix,2, median)

colnames(SCS)<-c("lower","upper")

out<-list(
conf.int=SCS,
estimate=estimate,
x=x,
k=k,
N=N,
conf.level=conf.level,
alternative=alternative)

class(out)<-"SCSnp"

return(out)

}


`SCSnp.CCRatio` <-
function(x,...)
{
args<-list(...)

args$x<-x$chains

out<-do.call("SCSnp.default", args)

out$x<-x

return(out)

}

`SCSnp.CCDiff` <-
function(x,...)
{
args<-list(...)

args$x<-x$chains

out<-do.call("SCSnp.default", args)

out$x<-x

return(out)

}

`SCSnp.bugs` <-
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

out<-do.call("SCSnp.default", args)

return(out)
}


