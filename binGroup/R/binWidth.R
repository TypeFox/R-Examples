"binWidth" <-
function(n, p, conf.level=0.95, alternative="two.sided", method="CP")
{

if( n < 3 ) {stop(" nmax must be at least 3 ")} 

if(length(n)<1 || length(n)>2 ) {stop(" n must be either a single integer or a vector of two integers")} 

if(conf.level<0 || conf.level>1 ) {stop(" conf.level must be a value between 0 and 1 ")} 

if(p<0 || p>1 ) {stop(" p.hyp must be a value between 0 and 1")} 
 
method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

# indicator function for the CI length at a special event
# in one sided case: length is defined as absolute difference between estimator and confidence bound

 L.Ind.bin<-function(y, n, p, conf.level, alternative, method)

 {
  switch(method,
  "Wald"={int<-binWald(y=y,n=n, conf.level=conf.level, alternative=alternative)},
  "Score"={int<-binWilson(y=y,n=n, conf.level=conf.level, alternative=alternative)},
  "AC"={int<-binAC(y=y,n=n, conf.level=conf.level, alternative=alternative)},
  "SOC"={int<-binSOC(y=y,n=n, conf.level=conf.level, alternative=alternative)},
  "CP"={int<-binCP(y=y,n=n, conf.level=conf.level, alternative=alternative)},
  "Blaker"={int<-binBlaker(y=y,n=n, conf.level=conf.level)}
  )

  if(alternative=="less")
    {CIlength <- int[[2]]-p}

  if(alternative=="greater")
    {CIlength <- p-int[[1]]}

  if(alternative=="two.sided")
    {CIlength <- int[[2]]-int[[1]]}
  CIlength
 }

# Probability of a single event, the binomial density:

 bin.prob <- function(y,n,p)
  {
   exp( lchoose(n,y) + y*log(p) + (n-y)*log(1-p) )
  }

#  calculate this for all possible events: 

yvec<-0:n

   Lvec<-numeric(length=length(yvec))   
   probvec<-numeric(length=length(yvec))

   for(i in 1:length(yvec))
    {Lvec[i]<-L.Ind.bin(y=yvec[i], n=n, p=p, conf.level=conf.level, alternative=alternative, method=method)
     probvec[i]<-bin.prob(y=yvec[i], n=n, p=p)
    }
  expCILength=sum(Lvec * probvec)

# E(X)= sum(Xi * prob(Xi))

out<-list(expCIWidth=expCILength, alternative=alternative, p=p, n=n)

class(out)<-"binWidth"
return(out)
}

