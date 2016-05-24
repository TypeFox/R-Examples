powerbinom <-
function(p, n, alpha=0.05, type="Dunnett", cmat=NULL, rhs=0, alternative=c("two.sided", "less", "greater"), ptype=c("global", "anypair", "allpair"), method="Wald", crit=NULL, ...)
{

if(any(c(n*p, n*(1-p))<5)){warning("For small n and extreme proportions this function might give misleading approximations.")}

if(length(alpha)!=1 || !is.numeric(alpha) | alpha<=0 | alpha>0.5){stop("alpha must be a single numeric value between 0 and 0.5")}
 ptype <- match.arg(ptype)
 alternative <- match.arg(alternative)
if( length(n)<2 || !(is.numeric(n)|is.integer(n)) ){stop("n must be a vector of sample sizes (at least length 2, containing numeric or integer values)")}
if( length(p)<2 || !is.numeric(p) || (any(p<0)| any(p>1)) ){stop("p must be a vector of proportions (at least length 2, containing numeric values between 0 and 1)")}
ngroup<-length(p)
if(length(n)!=ngroup){stop("length of vector n must equal the length of vector p")}

if(is.null(cmat)){cmat<-contrMat(n=n, type=type)}else{
if(!is.matrix(cmat) | !is.numeric(cmat)){stop("cmat must be a matrix with numeric entries")}
if(ncol(cmat)!=ngroup){stop("number of columns in cmat must equal the length of p")}
#        if(nrow(cmat)<1){stop("number of rows (i.e. number of comparisons) in cmat should be at least 2")}
}

M <- nrow(cmat)

if(!is.numeric(rhs) & !is.integer(rhs)){stop("rhs must be a (vector of) numeric value(s).")}
if(any(rhs >= 1) | any(rhs <= (-1))){stop("rhs is bounded between -1 and 1 for differences of proportions.")}
if(alternative=="two.sided" & any(rhs != 0)){warning("Power of shifted tests (i.e. rhs!=0) makes sense for one-sided alternatives, but alternative='two.sided' has been specified.")}
RHS<-rep(rhs, length.out=M)

switch(method,
"Wald"={
  Ltrue <- as.numeric(cmat%*%p)
  LExpectH1 <- cmat%*%p
  VarpH1 <- p*(1-p)/n
  LVarH1 <- (cmat^2) %*% (VarpH1)
  ExpTH1 <-(LExpectH1-RHS)/sqrt(LVarH1)
 },
"ADD1"={
  Ltrue <- as.numeric(cmat%*%p)
  PExpH1<-(p*n+0.5)/(n+1)
  LExpectH1 <- cmat%*%PExpH1
  VarpH1 <- PExpH1*(1-PExpH1)/(n+1)
  LVarH1 <- (cmat^2) %*% (VarpH1)
  ExpTH1 <- (LExpectH1-RHS)/sqrt(LVarH1)
 },
"ADD2"={
  Ltrue <- as.numeric(cmat%*%p)
  PExpH1<-(p*n+1)/(n+2)
  LExpectH1 <- cmat%*%PExpH1
  VarpH1 <- PExpH1*(1-PExpH1)/(n+2)
  LVarH1 <- (cmat^2) %*% (VarpH1)
  ExpTH1 <- (LExpectH1-RHS)/sqrt(LVarH1)
 })

# calculate the correlation matrix

CORR <- corrMatgen(CM=cmat, varp=VarpH1)

if(is.null(crit)){
switch(alternative,
two.sided={crit <- qmvnorm(p=1-alpha, tail =  "both.tails", sigma = CORR, ...)[["quantile"]]},
less={crit <- qmvnorm(p=1-alpha, tail =  "upper",  sigma = CORR, ...)[["quantile"]]},
greater={crit <- qmvnorm(p=1-alpha, tail =  "lower", sigma = CORR,...)[["quantile"]]})

}else{
if(!is.numeric(crit)| length(numeric)!=1){stop("crit must be a single numeric value")}
switch(alternative,
two.sided={if(crit<0){crit <- abs(crit); warning("Specified critical value 'crit' is negative, its absolute value will be used instead.")}},
less={if(crit>0){crit <- (-1)*crit; warning("Specified critical value 'crit' is positive, its negative value will be used instead.")}},
greater={if(crit<0){crit <- abs(crit); warning("Specified critical value 'crit' is negative, its absolute value will be used instead.")}})
}

ExpTeststat<-as.numeric(ExpTH1)
M<-length(ExpTeststat)

switch(EXPR=ptype,
"global"={

switch(EXPR=alternative,
"two.sided"={beta <- pmvnorm(lower=rep(-crit,M), upper=rep(crit,M), mean=ExpTeststat, sigma=CORR,...); whichHA <- which(abs(Ltrue-RHS) > 10*.Machine$double.eps)},
"less"={beta <- pmvnorm(lower=rep(crit,M), upper=rep(Inf, M), mean=ExpTeststat, sigma=CORR,...); whichHA <- which(Ltrue-RHS < (-10)*.Machine$double.eps)},
"greater"={beta <- pmvnorm(lower=rep(-Inf,M), upper=rep(crit,M), mean=ExpTeststat, sigma=CORR,...); whichHA <- which(Ltrue-RHS > 10*.Machine$double.eps)})
},

"anypair"={

switch(EXPR=alternative,
"two.sided"={
whichHA <- which(abs(Ltrue-RHS) > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvnorm(lower=rep(-crit, MHA), upper=rep(crit, MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA,whichHA],...)}
},
"less"={
whichHA <- which(Ltrue-RHS < (-10)*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvnorm(lower=rep(crit, MHA), upper=rep(Inf, MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA,whichHA],...)}
},
"greater"={
whichHA <- which(Ltrue-RHS > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvnorm(lower=rep(-Inf,MHA), upper=rep(crit, MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA,whichHA],...)}
})
},

"allpair"={

switch(EXPR=alternative,
"two.sided"={
whichHA <- which(abs(Ltrue-RHS) > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
nsim <- 10000
RT <- rmvnorm(n=nsim, mean=ExpTeststat[whichHA], sigma=CORR[whichHA, whichHA], method="svd")
nreject <- sum( apply(RT, 1, function(x){min(abs(x))}) > abs(crit) )
beta <- 1 - (nreject/nsim)
simerror <- sqrt(beta*(1-beta)/nsim)
attr(beta, which="simerror")<-simerror
}},
"less"={
whichHA <- which(Ltrue-RHS < (-10)*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
beta <- 1 - pmvnorm(lower=rep(-Inf, MHA), upper=rep(crit,MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA, whichHA],...)
}},
"greater"={
whichHA <- which(Ltrue-RHS > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
beta <- 1 - pmvnorm(lower=rep(crit,MHA), upper=rep(Inf,MHA), mean=ExpTeststat[whichHA],  sigma=CORR[whichHA, whichHA],...)
}})
})

pow<- (1-beta)

HAtrue<-numeric(length=M)
HAtrue[whichHA]<-1


if(is.null(colnames(cmat))){colnames(cmat)<-paste("g", 1:ngroup, sep="")}
settings<-data.frame(cmat, trueContrast=Ltrue, rhs=RHS, ExpTstat=ExpTeststat, underHA=HAtrue)

RES<-list(power=pow,
n=n,
p=p,
conexp=settings,
crit=crit,
alternative=alternative,
ptype=ptype,
alpha=alpha
)

return(RES)
}



powerbinomOR <- function (p, n, alpha = 0.05, type = "Dunnett", cmat = NULL, 
    rhs = 1, alternative = c("two.sided", "less", "greater"), 
    ptype = c("global", "anypair", "allpair"), 
    crit = NULL, ...) 
{
    if (any(c(n * p, n * (1 - p)) < 5)) {
        warning("For small n and extreme proportions this function might give misleading approximations.")
    }
    if (length(alpha) != 1 || !is.numeric(alpha) | alpha <= 0 | 
        alpha > 0.5) {
        stop("alpha must be a single numeric value between 0 and 0.5")
    }
    ptype <- match.arg(ptype)
    alternative <- match.arg(alternative)
    if (length(n) < 2 || !(is.numeric(n) | is.integer(n))) {
        stop("n must be a vector of sample sizes (at least length 2, containing numeric or integer values)")
    }
    if (length(p) < 2 || !is.numeric(p) || (any(p <= 0) | any(p >=
        1))) {
        stop("p must be a vector of proportions (at least length 2, containing numeric values greater than 0 and smaller than 1)")
    }
    ngroup <- length(p)
    if (length(n) != ngroup) {
        stop("length of vector n must equal the length of vector p")
    }
    if (is.null(cmat)) {
        cmat <- contrMat(n = n, type = type)
    }
    else {
        if (!is.matrix(cmat) | !is.numeric(cmat)) {
            stop("cmat must be a matrix with numeric entries")
        }
        if (ncol(cmat) != ngroup) {
            stop("number of columns in cmat must equal the length of p")
        }
    }
    M <- nrow(cmat)
    if (!is.numeric(rhs) & !is.integer(rhs)) {
        stop("rhs must be a (vector of) numeric value(s).")
    }
    if (any(rhs <= 0)) {
        stop("rhs must be greater 0 for odds ratios.")
    }
    if (alternative == "two.sided" & any(rhs != 1)) {
        warning("Power of shifted tests (i.e. rhs!=1) makes sense for one-sided alternatives, but alternative='two.sided' has been specified.")
    }
    RHS <- rep(rhs, length.out = M)
	LRHS <- log(RHS)
	logodds <- log(p/(1-p))
        Ltrue <- as.numeric(cmat %*% logodds)
        LExpectH1 <- cmat %*% logodds
        VarpH1 <- 1/(p*n) + 1/((1-p)*n)
        LVarH1 <- (cmat^2) %*% (VarpH1)
        ExpTH1 <- (LExpectH1 - LRHS)/sqrt(LVarH1)
 
    CORR <- corrMatgen(CM = cmat, varp = VarpH1)
    if (is.null(crit)) {
        switch(alternative, two.sided = {
            crit <- qmvnorm(p = 1 - alpha, tail = "both.tails", 
                sigma = CORR, ...)[["quantile"]]
        }, less = {
            crit <- qmvnorm(p = 1 - alpha, tail = "upper", sigma = CORR, 
                ...)[["quantile"]]
        }, greater = {
            crit <- qmvnorm(p = 1 - alpha, tail = "lower", sigma = CORR, 
                ...)[["quantile"]]
        })
    }
    else {
        if (!is.numeric(crit) | length(numeric) != 1) {
            stop("crit must be a single numeric value")
        }
        switch(alternative, two.sided = {
            if (crit < 0) {
                crit <- abs(crit)
                warning("Specified critical value 'crit' is negative, its absolute value will be used instead.")
            }
        }, less = {
            if (crit > 0) {
                crit <- (-1) * crit
                warning("Specified critical value 'crit' is positive, its negative value will be used instead.")
            }
        }, greater = {
            if (crit < 0) {
                crit <- abs(crit)
                warning("Specified critical value 'crit' is negative, its absolute value will be used instead.")
            }
        })
    }
    ExpTeststat <- as.numeric(ExpTH1)
    M <- length(ExpTeststat)
    switch(EXPR = ptype, global = {
        switch(EXPR = alternative, two.sided = {
            beta <- pmvnorm(lower = rep(-crit, M), upper = rep(crit, 
                M), mean = ExpTeststat, sigma = CORR, ...)
            whichHA <- which(abs(Ltrue - LRHS) > 10 * .Machine$double.eps)
        }, less = {
            beta <- pmvnorm(lower = rep(crit, M), upper = rep(Inf, 
                M), mean = ExpTeststat, sigma = CORR, ...)
            whichHA <- which(Ltrue < LRHS)
        }, greater = {
            beta <- pmvnorm(lower = rep(-Inf, M), upper = rep(crit, 
                M), mean = ExpTeststat, sigma = CORR, ...)
            whichHA <- which(Ltrue > LRHS)
        })
    }, anypair = {
        switch(EXPR = alternative, two.sided = {
            whichHA <- which(abs(Ltrue - LRHS) > 10 * .Machine$double.eps)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- pmvnorm(lower = rep(-crit, MHA), upper = rep(crit, 
                  MHA), mean = ExpTeststat[whichHA], sigma = CORR[whichHA, 
                  whichHA], ...)
            }
        }, less = {
            whichHA <- which(Ltrue < LRHS)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- pmvnorm(lower = rep(crit, MHA), upper = rep(Inf, 
                  MHA), mean = ExpTeststat[whichHA], sigma = CORR[whichHA, 
                  whichHA], ...)
            }
        }, greater = {
            whichHA <- which(Ltrue > LRHS)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- pmvnorm(lower = rep(-Inf, MHA), upper = rep(crit, 
                  MHA), mean = ExpTeststat[whichHA], sigma = CORR[whichHA, 
                  whichHA], ...)
            }
        })
    }, allpair = {
        switch(EXPR = alternative, two.sided = {
            whichHA <- which(abs(Ltrue - LRHS) > 10 * .Machine$double.eps)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                nsim <- 10000
                RT <- rmvnorm(n = nsim, mean = ExpTeststat[whichHA], 
                  sigma = CORR[whichHA, whichHA], method = "svd")
                nreject <- sum(apply(RT, 1, function(x) {
                  min(abs(x))
                }) > abs(crit))
                beta <- 1 - (nreject/nsim)
                simerror <- sqrt(beta * (1 - beta)/nsim)
                attr(beta, which = "simerror") <- simerror
            }
        }, less = {
            whichHA <- which(Ltrue < LRHS)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- 1 - pmvnorm(lower = rep(-Inf, MHA), upper = rep(crit, 
                  MHA), mean = ExpTeststat[whichHA], sigma = CORR[whichHA, 
                  whichHA], ...)
            }
        }, greater = {
            whichHA <- which(Ltrue > LRHS)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- 1 - pmvnorm(lower = rep(crit, MHA), upper = rep(Inf, 
                  MHA), mean = ExpTeststat[whichHA], sigma = CORR[whichHA, 
                  whichHA], ...)
            }
        })
    })
    pow <- (1 - beta)
    HAtrue <- numeric(length = M)
    HAtrue[whichHA] <- 1
    if (is.null(colnames(cmat))) {
        colnames(cmat) <- paste("g", 1:ngroup, sep = "")
    }
    settings <- data.frame(cmat, trueOddsRatio = exp(Ltrue), rhs = RHS, 
        trueLogit=Ltrue, rhsLogitscale=LRHS,
        ExpTstat = ExpTeststat, underHA = HAtrue)
    RES <- list(power = pow, n = n, p = p, conexp = settings, 
        crit = crit, alternative = alternative, ptype = ptype, 
        alpha = alpha)
    return(RES)
}


