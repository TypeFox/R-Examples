powermcpt <-
function(mu, n, sd, cmat=NULL, rhs=0, type="Dunnett", alternative=c("two.sided", "less", "greater"), alpha=0.05, ptype=c("global", "anypair", "allpair"), crit=NULL, ...)
{

rmvtFS <- function (n, sigma = diag(2), df = 1, delta = rep(0, nrow(sigma)), 
    type = c("shifted", "Kshirsagar"), method=c("eigen", "svd", "chol")) 
{
    if (length(delta) != nrow(sigma)) 
        stop("delta and sigma have non-conforming size")
    if (df == 0) 
        return(rmvnorm(n, mean = delta, sigma = sigma, method=method))
    type <- match.arg(type)
    if (type == "Kshirsagar") 
        return(rmvnorm(n, mean = delta, sigma = sigma, method=method)/sqrt(rchisq(n, 
            df)/df))
    if (type == "shifted") {
        sims <- rmvnorm(n, sigma = sigma, method=method)/sqrt(rchisq(n, df)/df)
        return(sweep(sims, 2, delta, "+"))
    }
}

if(length(sd)!=1 || !is.numeric(sd)){stop("sd must be a single numeric value")}
if(length(alpha)!=1 || !is.numeric(alpha) | alpha<=0 | alpha>0.5){stop("alpha must be a single numeric value between 0 and 0.5")}
 ptype <- match.arg(ptype)
 alternative <- match.arg(alternative)
if( length(n)<2 || !(is.numeric(n)|is.integer(n)) ){stop("n must be a vector of sample sizes (at least length 2, containing numeric or integer values)")}
if( length(mu)<2 || !(is.numeric(mu)|is.integer(mu)) ){stop("mu must be a vector of expected means (at least length 2, containing numeric or integer values)")}
        ngroup<-length(mu)
        if(length(n)!=ngroup){stop("length of vector n must equal the length of vector mu")}

if(is.null(cmat)){cmat<-contrMat(n=n, type=type)}else{
if(!is.matrix(cmat) | !is.numeric(cmat))stop("cmat must be a matrix with numeric entries")
        if(ncol(cmat)!=ngroup){stop("number of columns in cmat must equal the length of mu")}
        if(nrow(cmat)<1){stop("number of rows (i.e. number of comparisons) in cmat should be at least 2")}
}

if(!is.numeric(rhs) & !is.integer(rhs)){stop("rhs must be a (vector of) numeric value(s).")}

M <- nrow(cmat)
MU <- matrix(mu, ncol=1)
VCOV <- diag((sd^2)/n)
VCOVL <- cmat %*% VCOV %*% t(cmat)
CORR <- cov2cor(VCOVL)
dfR <- sum(n)-length(n)

RHS<-rep(rhs, length.out=M)


if(is.null(crit)){
switch(alternative,
two.sided={crit <- qmvt(p=1-alpha,tail =  "both.tails", df = dfR, corr = CORR, ...)[["quantile"]]},
less={crit <- qmvt(p=1-alpha,tail =  "upper", df = dfR,  corr = CORR, ...)[["quantile"]]},
greater={crit <- qmvt(p=1-alpha,tail =  "lower", df = dfR, corr = CORR,...)[["quantile"]]})

}else{
if(!is.numeric(crit)| length(numeric)!=1){stop("crit must be a single numeric value")}
switch(alternative,
two.sided={if(crit<0){crit <- abs(crit); warning("Specified critical value 'crit' is negative, its absolute value will be used instead.")}},
less={if(crit>0){crit <- (-1)*crit; warning("Specified critical value 'crit' is positive, its negative value will be used instead.")}},
greater={if(crit<0){crit <- abs(crit); warning("Specified critical value 'crit' is negative, its absolute value will be used instead.")}})
}


ExpL <- cmat %*% MU
Ltrue <- as.numeric(ExpL)
ExpVL <- diag(VCOVL)
ExpTeststat <- (as.numeric(ExpL)-RHS)/sqrt(as.numeric(ExpVL))

switch(EXPR=ptype,
"global"={

switch(EXPR=alternative,
"two.sided"={beta <- pmvt(lower=rep(-crit,M), upper=rep(crit,M), delta=ExpTeststat, df=dfR, corr=CORR,...); whichHA <- which(abs(Ltrue-RHS) > 10*.Machine$double.eps)},
"less"={beta <- pmvt(lower=rep(crit,M), upper=rep(Inf, M), delta=ExpTeststat, df=dfR, corr=CORR,...); whichHA <- which(Ltrue-RHS < (-10)*.Machine$double.eps)},
"greater"={beta <- pmvt(lower=rep(-Inf,M), upper=rep(crit,M), delta=ExpTeststat, df=dfR, corr=CORR,...); whichHA <- which(Ltrue-RHS > 10*.Machine$double.eps)})
},

"anypair"={

switch(EXPR=alternative,
"two.sided"={
whichHA <- which(abs(Ltrue-RHS) > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvt(lower=rep(-crit, MHA), upper=rep(crit, MHA), delta=ExpTeststat[whichHA], df=dfR, corr=CORR[whichHA,whichHA],...)}
},
"less"={
whichHA <- which(Ltrue-RHS < (-10)*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvt(lower=rep(crit, MHA), upper=rep(Inf, MHA), delta=ExpTeststat[whichHA], df=dfR, corr=CORR[whichHA,whichHA],...)}
},
"greater"={
whichHA <- which(Ltrue-RHS > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvt(lower=rep(-Inf,MHA), upper=rep(crit, MHA), delta=ExpTeststat[whichHA], df=dfR, corr=CORR[whichHA,whichHA],...)}
})
},

"allpair"={

switch(EXPR=alternative,
"two.sided"={
whichHA <- which(abs(Ltrue-RHS) > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
nsim <- 10000
RT <- rmvtFS(n=nsim, delta=ExpTeststat[whichHA], df=dfR, sigma=CORR[whichHA, whichHA], method="svd")
nreject <- sum( apply(RT, 1, function(x){min(abs(x))}) > abs(crit) )
beta <- 1 - (nreject/nsim)
simerror <- sqrt(beta*(1-beta)/nsim)
attr(beta, which="simerror")<-simerror
}},
"less"={
whichHA <- which(Ltrue-RHS < (-10)*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
beta <- 1 - pmvt(lower=rep(-Inf, MHA), upper=rep(crit,MHA), delta=ExpTeststat[whichHA], df=dfR, corr=CORR[whichHA, whichHA],...)
}},
"greater"={
whichHA <- which(Ltrue-RHS > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
beta <- 1 - pmvt(lower=rep(crit,MHA), upper=rep(Inf,MHA), delta=ExpTeststat[whichHA], df=dfR, corr=CORR[whichHA, whichHA],...)
}})
})


pow<- (1-beta)

if(is.null(colnames(cmat))){colnames(cmat)<-paste("g", 1:ngroup, sep="")}

HAtrue<-numeric(length=M)
HAtrue[whichHA]<-1

settings<-data.frame(cmat, expContrast=ExpL, rhs=RHS, ExpTstat=ExpTeststat, underHA=HAtrue)

out<-list(power=pow,
mu=mu, n=n,
conexp=settings,
crit=crit,
alternative=alternative,
ptype=ptype,
alpha=alpha
)

return(out)
}

