"binAC" <-
function(n, y, conf.level=0.95, alternative="two.sided")
{
alpha=1-conf.level
est=y/n
z1s=qnorm(conf.level)
z2s=qnorm(1-alpha/2)

esti1s=(y+(z1s^2)/2)/(n+z1s^2)
esti2s=(y+(z2s^2)/2)/(n+z2s^2)

ni1s=n+z1s^2
ni2s=n+z2s^2

if(alternative=="two.sided"){

CI=c(esti2s-z2s*sqrt(esti2s*(1-esti2s)/(ni2s)),
     esti2s+z2s*sqrt(esti2s*(1-esti2s)/(ni2s)) )
}
else{if (alternative=="less"){
CI=c( 0 , esti1s+z1s*sqrt(esti1s*(1-esti1s)/(ni1s)) )
}

else{if(alternative=="greater"){
CI=c(esti1s-z1s*sqrt(esti1s*(1-esti1s)/(ni1s)), 1 )
}
else {stop("alternative mis-specified")}}}

CI
}

"binBlaker" <-
function (n,y,conf.level=0.95, tolerance=1e-04, alternative="two.sided")
{
acceptbin <- function(y,n,p)
{
  p1 = 1-pbinom(y-1, n, p)
  p2 = pbinom(y, n, p)
  a1 = p1 + pbinom( qbinom(p1,n,p)-1, n, p )
  a2 = p2+1-pbinom( qbinom(1-p2,n,p), n, p )
  return(min(a1,a2))
}

lower<-0
upper<-1

if(y!=0)
  {lower<-qbeta((1-conf.level)/2, y, n-y+1)
    {while(acceptbin(y,n,lower+tolerance)<(1-conf.level))
    lower=lower+tolerance}
  }

if(y!=n)
  {upper<-qbeta(1-(1-conf.level)/2, y+1, n-y)
    {while(acceptbin(y,n,upper-tolerance)<(1-conf.level))
    upper=upper-tolerance}
  }

c(lower, upper)

}

"binCP" <-
function(n, y, conf.level=0.95, alternative="two.sided")
{
lower<-0
upper<-1
if(alternative=="two.sided")
{
 if(y!=0)
   {lower<-qbeta((1-conf.level)/2, y, n-y+1)}
 if(y!=n)
   {upper<-qbeta(1-(1-conf.level)/2, y+1, n-y)}
}
if(alternative=="less")
{
 if(y!=n)
   {upper<-qbeta(1-(1-conf.level), y+1, n-y)}
}
if(alternative=="greater")
{
 if(y!=0)
   {lower<-qbeta((1-conf.level), y, n-y+1)}
}

c(lower,upper)

}


"binCI" <-
function(n, y, conf.level=0.95, alternative="two.sided", method="CP")

{
if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
if(length(y)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
if( length(conf.level)!=1 || conf.level<0 || conf.level>1){stop("conf.level must be a positive number between 0 and 1, usually 0.95")}

method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

estimate=y/n

switch(method,
"CP"={conf.int<-binCP(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"Blaker"={
 if(alternative=="less" || alternative=="greater"){
  warning("The Blaker CI is inherently two-sided.")
  conf.int<-c(NA, NA)
  }
 if(alternative=="two.sided")
  {conf.int<-binBlaker(n=n, y=y, conf.level=conf.level)}
 },
 
"Score"={conf.int<-binWilson(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"AC"={conf.int<-binAC(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"Wald"={conf.int<-binWald(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"SOC"={conf.int<-binSOC(n=n, y=y, conf.level=conf.level, alternative=alternative)}
)

out<-list(conf.int=conf.int,
	estimate=estimate,
	method=method,
	conf.level=conf.level,
	alternative=alternative)
class(out)<-"binCI"
out
}


"binSOC" <-
function(n,y,conf.level=0.95,alternative="two.sided")

{
esti<-y/n
kappa<-qnorm(conf.level)
eta<-(kappa^2)/3 + 1/6
gamma1<-((13/18)*kappa^2 + 17/18)*(-1)
gamma2<-(kappa^2)/18 + 7/36

midpo<-(y+eta)/(n+2*eta)

if(alternative=="less")
  {CI=c( 0 , midpo + kappa * sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n) )}

if(alternative=="greater")
  {CI=c( midpo - kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n) , 1)}

if (alternative=="two.sided")
{
kappa<-qnorm(1-(1-conf.level)/2)
eta<-(kappa^2)/3 + 1/6
gamma1<-((13/18)*kappa^2 + 17/18)*(-1)
gamma2<-(kappa^2)/18 + 7/36

CI=c( midpo - kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n) , 
 midpo + kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n) )
}
CI
}

"binWald" <-
function(n, y, conf.level=0.95, alternative="two.sided")
{
alpha=1-conf.level
est=y/n
z1s=qnorm(conf.level)
z2s=qnorm(1-alpha/2)

if(alternative=="two.sided"){
KI=c(est-z2s*sqrt(est*(1-est)/(n)),
     est+z2s*sqrt(est*(1-est)/(n)) )
}
else{if (alternative=="less"){
KI=c( 0 , est+z1s*sqrt(est*(1-est)/(n)) )
}
else{if(alternative=="greater"){
KI=c(est-z1s*sqrt(est*(1-est)/(n)), 1 )
}
else {stop("alternative mis-specified")}}}
KI
}

"binWilson" <-
function(n,y,conf.level=0.95,alternative="two.sided") {
alpha=1-conf.level
t=y/n
if(alternative =="two.sided"){
	est.int=(y+(qnorm(1-alpha/2)^2)/2)/(n+(qnorm(1-alpha/2))^2)
	w.se=((qnorm(1-alpha/2))*sqrt(n*t*(1-t)+(qnorm(1-alpha/2)^2)/4))/(n+qnorm(1-alpha/2)^2)
	CI=c( est.int-w.se, est.int+w.se )
}
else{if(alternative=="less"){
	est.int=(y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)
	w.se=((qnorm(1-alpha))*sqrt(n*t*(1-t)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
	CI=c( 0, est.int+w.se )
} 
else{if(alternative=="greater"){
	est.int=(y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)
	w.se=((qnorm(1-alpha))*sqrt(n*t*(1-t)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
	CI=c( est.int-w.se , 1 )
}
else{stop("argument alternative misspecified")}}}

CI

}


