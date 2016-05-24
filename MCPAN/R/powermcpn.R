powermcpn <-
function(ExpTeststat, corrH1, crit, alternative=c("two.sided", "less", "greater"), alpha=0.05, ptype=c("global", "anypair", "allpair"), ...)
{

M<-length(ExpTeststat)
CORR<-corrH1

switch(EXPR=ptype,
"global"={

switch(EXPR=alternative,
"two.sided"={beta <- pmvnorm(lower=rep(-crit,M), upper=rep(crit,M), mean=ExpTeststat, sigma=CORR,...); whichHA <- which(abs(ExpTeststat) > 10*.Machine$double.eps)},
"less"={beta <- pmvnorm(lower=rep(crit,M), upper=rep(Inf, M), mean=ExpTeststat, sigma=CORR,...); whichHA <- which(ExpTeststat < -10*.Machine$double.eps)},
"greater"={beta <- pmvnorm(lower=rep(-Inf,M), upper=rep(crit,M), mean=ExpTeststat, sigma=CORR,...); whichHA <- which(ExpTeststat > 10*.Machine$double.eps)})
},

"anypair"={

switch(EXPR=alternative,
"two.sided"={
whichHA <- which(abs(ExpTeststat) > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvnorm(lower=rep(-crit, MHA), upper=rep(crit, MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA,whichHA],...)}
},
"less"={
whichHA <- which(ExpTeststat < -10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvnorm(lower=rep(crit, MHA), upper=rep(Inf, MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA,whichHA],...)}
},
"greater"={
whichHA <- which(ExpTeststat > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated."); beta <- 1-alpha}else{
beta <- pmvnorm(lower=rep(-Inf,MHA), upper=rep(crit, MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA,whichHA],...)}
})
},

"allpair"={

switch(EXPR=alternative,
"two.sided"={
whichHA <- which(abs(ExpTeststat) > 10*.Machine$double.eps)
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
whichHA <- which(ExpTeststat < -10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
beta <- 1 - pmvnorm(lower=rep(-Inf, MHA), upper=rep(crit,MHA), mean=ExpTeststat[whichHA], sigma=CORR[whichHA, whichHA],...)
}},
"greater"={
whichHA <- which(ExpTeststat > 10*.Machine$double.eps)
MHA <- length(whichHA)
if(MHA<1){warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated."); beta<-1-alpha}else{
beta <- 1 - pmvnorm(lower=rep(crit,MHA), upper=rep(Inf,MHA), mean=ExpTeststat[whichHA],  sigma=CORR[whichHA, whichHA],...)
}})
})


pow<- (1-beta)

HAtrue<-numeric(length=M)
HAtrue[whichHA]<-1

settings<-data.frame(ExpTstat=ExpTeststat, underHA=HAtrue)

out<-list(power=pow,
conexp=settings,
crit=crit,
alternative=alternative,
ptype=ptype,
alpha=alpha
)

return(out)
}

