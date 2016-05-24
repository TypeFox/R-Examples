"binTest" <-
function(n, y, p.hyp, alternative="two.sided", method="Exact")
{

if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
if(length(y)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}

if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
if(length(p.hyp)!=1 || p.hyp<0 || p.hyp>1){stop("p.hyp must be a positive number between 0 and 1")}

method<-match.arg(method, choices=c("Exact","Score","Wald"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

 esti=y/n 

switch(method,
"Exact"={p.val=binom.test(x=y,n=n,p=p.hyp, alternative=alternative)$p.value},
 
"Score"={
   teststat = (esti-p.hyp)/(sqrt(p.hyp*(1-p.hyp)/n))
   if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
   if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
   if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
  },

"Wald"={
   teststat = (esti-p.hyp)/(sqrt(esti*(1-esti)/n))
   if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
   if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
   if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)}  
  }
)
 
out<-list(p.value=p.val,
estimate=esti,
alternative=alternative,
p.hyp=p.hyp,
method=method)

class(out)<-"binTest"
out
}

