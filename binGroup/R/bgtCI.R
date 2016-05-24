"bgtAC" <-
function(n, y, s, conf.level=0.95, alternative="two.sided") {
alpha=1-conf.level

est.int=(y+(qnorm(1-alpha/2)^2)/2)/(n+(qnorm(1-alpha/2))^2)
est.int1s=(y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)


if(alternative =="two.sided"){
AC.se=(qnorm(1-alpha/2))*sqrt((est.int*(1-est.int))/(n+(qnorm(1-alpha/2))^2))
KI.int.l=est.int-AC.se
KI.int.u=est.int+AC.se
        if (KI.int.u>1){KI.int.u=1}
        if (KI.int.l<0){KI.int.l=0}
CI=c(1-(1-KI.int.l)^(1/s),1-(1-KI.int.u)^(1/s))
} 

else{if(alternative=="less"){
AC.se=(qnorm(1-alpha))*sqrt((est.int1s*(1-est.int1s))/(n+(qnorm(1-alpha))^2))
KI.int.u=est.int1s+AC.se
        if (KI.int.u>1){KI.int.u=1}
CI=c(0, 1-(1-KI.int.u)^(1/s))
} 

else{if(alternative=="greater"){
AC.se=(qnorm(1-alpha))*sqrt((est.int1s*(1-est.int1s))/(n+(qnorm(1-alpha))^2))
KI.int.l=est.int1s-AC.se
        if (KI.int.l<0){KI.int.l=0}
CI=c(1-(1-KI.int.l)^(1/s),1)
}

else{stop("argument alternative misspecified")}}}

CI
}

"bgtBlaker" <-
function (n, y, s, conf.level=0.95, alternative="two.sided")
{

# # # from the S code given in Blaker(2000), slightly changed # # #

tolerance=1e-04

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
CI=c( 1-(1-lower)^(1/s), 1-(1-upper)^(1/s) )

CI  
}



"bgtCI" <-
function(n, s, y, conf.level=0.95, alternative="two.sided", method="CP")

{
if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
if(length(s)!=1 || (s<1 | abs(round(s)-s) > 1e-07)){stop("group size s must be specified as a single integer > 0")}
if(length(s)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
if(length(conf.level)!=1 || conf.level<0 || conf.level>1){stop("conf.level must be a positive number between 0 and 1")}

method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

estimate=1-(1-y/n)^(1/s)

switch(method,
"CP"={conf.int<-bgtCP(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"Blaker"={
 if(alternative=="less" || alternative=="greater"){
    warning("The Blaker CI is inherently two.sided")
    conf.int<-c(NA, NA)}
 if(alternative=="two.sided")
   {conf.int<-bgtBlaker(n=n, s=s, y=y, conf.level=conf.level)}
 },

"Score"={conf.int<-bgtWilson(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"AC"={conf.int<-bgtAC(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"Wald"={conf.int<-bgtWald(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"SOC"={conf.int<-bgtSOC(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)}
)

out<-list(conf.int=conf.int,
	estimate=estimate,
	method=method,
	conf.level=conf.level,
	alternative=alternative)
class(out)<-"bgtCI"
out
}

"bgtCP" <-
function(n, y, s, conf.level=0.95, alternative="two.sided")
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

estimate=1-(1-y/n)^(1/s)


CI=c(1-(1-lower)^(1/s),1-(1-upper)^(1/s))

CI   
}

"bgtSOC" <-
function(n,s,y,conf.level=0.95,alternative="two.sided")

{
esti<-y/n
kappa<-qnorm(conf.level)
eta<-(kappa^2)/3 + 1/6
gamma1<-((13/18)*kappa^2 + 17/18)*(-1)
gamma2<-(kappa^2)/18 + 7/36

midpo<-(y+eta)/(n+2*eta)

if(alternative=="less")
  {upper = midpo + kappa * sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n)
  CI=c( 0 ,upper)
  if(y==n||upper>1){CI=c(0,1)}
  else{ CI=c( 0 ,upper)}
 }

if(alternative=="greater")
  {CI=c( midpo - kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n) , 1)
  if(y==0){CI=c(0,1)} }

if (alternative=="two.sided")
{
kappa<-qnorm(1-(1-conf.level)/2)
eta<-(kappa^2)/3 + 1/6
gamma1<-((13/18)*kappa^2 + 17/18)*(-1)
gamma2<-(kappa^2)/18 + 7/36

lower= midpo - kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n)  
upper= midpo + kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n)
 
if(y==0){CI=c(0,upper)} 
else{if(y==n||upper>1){CI=c(lower,1)}
else{CI=c(lower, upper)}}
}

CI2=c( 1-(1-CI[1])^(1/s) ,  1-(1-CI[2])^(1/s))
CI2

}

"bgtWald" <-
function (n, y, s, conf.level=0.95, alternative="two.sided")

{
if(y>n) {stop("number of positive tests y can not be greater than number of tests n")}
th=y/n
esti=1-(1-th)^(1/s)
var.esti=(1-(1-esti)^s)/(n*(s^2)*(1-esti)^(s-2))
alpha=1-conf.level

if(alternative=="two.sided"){
    snquant=qnorm(p=1-alpha/2,mean=0,sd=1,lower.tail=TRUE)
    CI=c(esti-snquant*sqrt(var.esti),esti+snquant*sqrt(var.esti))
}
else{if (alternative=="less"){
    snquant=qnorm(p=1-alpha,mean=0,sd=1,lower.tail=TRUE)
    CI=c(0 ,esti+snquant*sqrt(var.esti))
}
else {if (alternative=="greater"){
    snquant=qnorm(p=1-alpha,mean=0,sd=1,lower.tail=TRUE)
    CI=c(esti-snquant*sqrt(var.esti), 1)
}
else {stop("argument alternative mis-specified")}}}
CI
}

"bgtWilson" <-
function(n, y, s, conf.level=0.95, alternative="two.sided")
{ 
alpha=1-conf.level 
th=y/n
est.int=(y+(qnorm(1-alpha/2)^2)/2)/(n+(qnorm(1-alpha/2))^2)
est.int1s=(y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)

if(alternative =="two.sided"){
    w.se=((qnorm(1-alpha/2))*sqrt(n*th*(1-th)+(qnorm(1-alpha/2)^2)/4))/(n+qnorm(1-alpha/2)^2)
    KI.int.l=est.int-w.se
    KI.int.u=est.int+w.se
        if (KI.int.u>1){KI.int.u=1}
        if (KI.int.l<0){KI.int.l=0}
    KI=c( 1-(1-KI.int.l)^(1/s), 1-(1-KI.int.u)^(1/s) )
}

else{if(alternative=="less"){
    w.se=((qnorm(1-alpha))*sqrt(n*th*(1-th)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
    KI.int.u=est.int1s+w.se
        if (KI.int.u>1){KI.int.u=1}
    KI=c( 0, 1-(1-KI.int.u)^(1/s) )
}

else{if(alternative=="greater"){
    w.se=((qnorm(1-alpha))*sqrt(n*th*(1-th)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
    KI.int.l=est.int1s-w.se
        if (KI.int.l<0){KI.int.l=0}
    KI=c( 1-(1-KI.int.l)^(1/s), 1 )
}

else{stop("argument alternative misspecified")}}}

KI

}









