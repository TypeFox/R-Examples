## ------------------------------------------------------------------------
require(softImpute)
set.seed(1011)
x=matrix(rnorm(30),6,5)
x[sample(1:30,10,replace=FALSE)]=NA
x
fits=softImpute(x,trace=TRUE,type="svd")
fits

## ------------------------------------------------------------------------
fita=softImpute(x,trace=TRUE)

## ------------------------------------------------------------------------
fits2=softImpute(x,rank.max=3,lambda=1.9,trace=TRUE,type="svd")
fita2=softImpute(x,rank.max=3,lambda=1.9,trace=TRUE)
fits2$d

## ------------------------------------------------------------------------
complete(x,fits2)

## ------------------------------------------------------------------------
xc=biScale(x,col.scale=FALSE,row.scale=FALSE,trace=TRUE)
xc
fits3=softImpute(xc,rank.max=3,lambda=1,type="svd")
fits3$d
complete(x,fits3,unscale=TRUE)

## ------------------------------------------------------------------------
fits2$d
deBias(x,fits2)$d

## ------------------------------------------------------------------------
xs=as(x,"Incomplete")
xs

## ------------------------------------------------------------------------
i=row(x)[!is.na(x)]
j=col(x)[!is.na(x)]
value=x[!is.na(x)]
cbind(i,j,value)

## ------------------------------------------------------------------------
Incomplete(i=i,j=j,x=value)

## ------------------------------------------------------------------------
xsc=biScale(xs,col.scale=FALSE,row.scale=FALSE)
fitss=softImpute(xsc,rank.max=3,lambda=1,trace=TRUE,type="svd")
fitss$d
fits3$d

## ------------------------------------------------------------------------
impute(fitss,i=c(2,3),j=c(2,2))

## ------------------------------------------------------------------------
x0=sparseMatrix(i=i,j=j,x=value)
x0
x0c=biScale(x0,col.scale=FALSE,row.scale=FALSE,row.center=FALSE)
x0c

## ------------------------------------------------------------------------
svdx0c=svd.als(x0c,rank=3,trace=TRUE)
svdx0c$d

## ------------------------------------------------------------------------
x02=as.matrix(x0)
svd(scale(x02,TRUE,FALSE))$d

## ------------------------------------------------------------------------
lam0=lambda0(xs)
lam0
fit0=softImpute(xs,lambda=lam0+.2)
fit0$d

## ------------------------------------------------------------------------
xs0=as(xs,"sparseMatrix")
fit0=svd.als(xs0)
fit0$d

## ------------------------------------------------------------------------
lamseq=exp(seq(from=log(lam0),to=log(1),length=10))
lamseq

## ------------------------------------------------------------------------
fits=as.list(lamseq)
ranks=as.integer(lamseq)
rank.max=2
warm=NULL
for( i in seq(along=lamseq)){
  fiti=softImpute(xs,lambda=lamseq[i],rank=rank.max,warm=warm)
  ranks[i]=sum(round(fiti$d,4)>0)
  rank.max=min(ranks[i]+2,4)
  warm=fiti
  fits[[i]]=fiti
  cat(i,"lambda=",lamseq[i],"rank.max",rank.max,"rank",ranks[i],"\n")
  }

