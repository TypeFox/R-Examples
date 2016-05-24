"bgtTest" <-
function(n, y, s, p.hyp, alternative="two.sided", method="Exact")
{

if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
if(length(s)!=1 || (s<1 | abs(round(s)-s) > 1e-07)){stop("group size s must be specified as a single integer > 0")}
if(length(s)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
if(length(p.hyp)!=1 || p.hyp<0 || p.hyp>1){stop("the proportion in the hypothesis p.hyp must be specified as a single value between 0 and 1")}

method<-match.arg(method, choices=c("Exact", "Score", "Wald"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

estimate=1-(1-y/n)^(1/s)

bgtsumProb<-function(x, n, s, p)
{
sumprob=0
for(i in x)
  {
  sumprob = sumprob + choose(n=n,k=i) * ((1-(1-p)^s)^(i)) * (1-p)^(s*(n-i))
  }
sumprob
}

switch(method,

"Exact"={
 if(alternative=="less")
  {p.val = bgtsumProb(x=0:y, n=n, s=s, p=p.hyp)}
 if(alternative=="greater")
  {p.val = bgtsumProb(x=y:n, n=n, s=s, p=p.hyp)}
 if (alternative=="two.sided")
  {p.val = min( 2*(bgtsumProb(x=0:y, n=n, s=s, p=p.hyp)), 2*(p.val = bgtsumProb(x=y:n, n=n, s=s, p=p.hyp)), 1)}
},


"Wald"={
esti = 1-(1-y/n)^(1/s)
varesti = (1-(1-esti)^s)/(n*(s^2)*(1-esti)^(s-2))  # variance estimator (see Swallow, 1985)
teststat = (esti-p.hyp)/sqrt(varesti) 

 if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
 if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
 if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
},

"Score"={
esti = y/n
t.hyp = 1-(1-p.hyp)^s
teststat = (esti-t.hyp)/(sqrt(t.hyp*(1-t.hyp)/n))

 if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
 if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
 if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
})

out<-list(p.value=p.val,
estimate=estimate,
alternative=alternative,
p.hyp=p.hyp,
method=method)

class(out)<-"bgtTest"
out
}

