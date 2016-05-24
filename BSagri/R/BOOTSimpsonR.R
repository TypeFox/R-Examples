`BOOTSimpsonR` <-
function(X, f, type="Dunnett",
 cmat=NULL, conf.level=0.95, alternative=c("two.sided", "less", "greater"), madj=TRUE, ...)
{
args<-list(...)
alternative<-match.arg(alternative)

BSimpson<-function(X, i, f)
{
XNEW<-as.data.frame(X[i,])
est<-estSimpsonf(X=XNEW, f=f)
return(est$estimate)
}

bargs<-args
bargs$data<-as.data.frame(X)
bargs$statistic=BSimpson
bargs$strata=f
bargs$f<-f
if(is.null(bargs$R)){bargs$R<-999}
if(is.null(bargs$sim)){bargs$sim<-"ordinary"}
if(is.null(bargs$stype)){bargs$stype<-"i"}

bootout<-do.call("boot", bargs)

ratiochains<-CCRatio.boot(x=bootout, cmat=cmat, type=type)

if(madj)
{
confint<-SCSnp.default(x=ratiochains$chains,
 conf.level=conf.level,
 alternative=alternative)
}else{
confint<-CInp.default(x=ratiochains$chains,
 conf.level=conf.level,
 alternative=alternative)
}

return(confint)

}

