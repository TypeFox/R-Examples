SCSrank <-
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

colnames(SCS)<-c("lower","upper")

attr(SCS, which="k")<-k
attr(SCS, which="N")<-N
OUT<-list(conf.int=SCS, conf.level=conf.level, alternative=alternative)
return(OUT)
}

