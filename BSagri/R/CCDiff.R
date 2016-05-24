`CCDiff` <-
function(bugs, dat, cmat=NULL,
 type=c("Dunnett", "Tukey", "Sequen", "Williams", "Changepoint"))
{

type<-match.arg(type)

if(class(bugs)!="bugs")
 {stop("argument bugs must be an object of class 'bugs'")}

if(class(dat)!="R2Bugsdat1w")
 {stop("argument dat must be an object of class 'R2Bugsdat1w'")}

if(dat$Intercept==TRUE)
 {stop("dat$Intercept must be FALSE")}

ngroup<-dat$names$ni

chains<-bugs$sims.list$muvec

if(is.null(cmat))
{
cmat<-contrMat(n=ngroup,type=type)
}
else{

if(!is.matrix(cmat))
 {stop("'cmat' must be a matrix, specifying the contrast coefficients")}

if(ngroup!=ncol(cmat))
 {stop("ncol(cmat) must be the same as the number of means in muvec")}

cs<-apply(cmat,1,sum)

if(any(cs!=0))
 {warning("Rows of cmat do not sum up to zero. Are the contrasts appropriately defined?")}

}

nchains<-apply(X=chains, MARGIN=1, FUN=function(x){cmat %*% x})

if(nrow(cmat)==1)
 {nchains<-matrix(nchains, nrow=1)}

rownames(nchains)<-rownames(cmat)

out<-list(
chains=t(nchains),
bugs=bugs,
dat=dat,
cmat=cmat
)

class(out)<-"CCDiff"
return(out)

}


`CCDiff.default` <-
function(x, cmat)
{


if(!is.matrix(x) & !is.data.frame(x))
 {stop("Argument 'x'must be a matrix or data.frame!")}

ngroup<-ncol(x)

Nsim<-nrow(x)

chains<-x

if(!is.matrix(cmat))
 {stop("'cmat' must be a matrix, specifying the contrast coefficients")}

if(ngroup!=ncol(cmat))
 {stop("ncol(cmat) must be the same as the number of means in muvec")}

cs<-apply(cmat,1,sum)

if(any(cs!=0))
 {warning("Rows of cmat do not sum up to zero. Are the contrasts appropriately defined?")}

nchains<-apply(X=chains, MARGIN=1, FUN=function(x){cmat %*% x})

if(nrow(cmat)==1)
 {nchains<-matrix(nchains, nrow=1)}

rownames(nchains)<-rownames(cmat)

out<-list(
chains=t(nchains),
x=x,
cmat=cmat
)

class(out)<-"CCDiff"
return(out)

}



`CCDiff.boot` <-
function(x, cmat=NULL,
 type=c("Dunnett","Tukey","Sequen","Williams","Changepoint","McDermott","GrandMean","Marcus"))
{

type<-match.arg(type)

if(type %in% c("Williams","Changepoint","McDermott","Marcus","GrandMean"))
 {warning("This is a test version. Choosing contrasts types differing from 'Dunnett','Tukey' or 'Sequen' might make no sense in case of unbalanced designs!")}

ngroup<-ncol(x$t)

f<-x$strata

ni<-unlist(lapply(split(f,f=f),length))

gnames<-names(x$t0)

names(ni)<-gnames

if(any(ni<5))
 {warning("For sample sizes les than 5 this function hardly makes sense!")}


if(is.null(cmat))
{
cmat<-contrMat(n=ni,type=type)
}

chains <- x$t

out<-CCDiff.default(x=chains, cmat=cmat)

return(out)

}




