"binPower" <-
function(n, delta, p.hyp, conf.level=0.95, alternative="two.sided", method="CP")
{

# Warnings:

if( n < 3 ) {stop(" nmax must be at least 3 ")} 

if(length(n)<1 || length(n)>2 ) {stop(" n must be either a single integer or a vector of two integers")} 

if(conf.level<0 || conf.level>1 ) {stop(" conf.level must be a value between 0 and 1 ")} 

if(p.hyp<0 || p.hyp>1 ) {stop(" p.hyp must be a value between 0 and 1")} 
 
if( delta<=0) {stop(" specify delta as absolute difference to p.hyp, thus as value greater than 0 ")} 

method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

 if(alternative=="less") 
  {
  if( p.hyp-delta <= 0 || p.hyp-delta >= 1 )
   {stop(" For alternative 'less': delta must be a positive value between 0 and p.hyp")} 
 
  }

 if(alternative=="greater") 
  {
  if( p.hyp+delta <= 0 || p.hyp+delta >= 1 ) 
   {stop("For alternative 'greater': delta must be a positive value between 0 and 1-p.hyp")} 
 
  }

 if(alternative=="two.sided") 

  {
  if( p.hyp+delta <= 0 || p.hyp+delta >= 1 || p.hyp-delta <= 0 || p.hyp-delta >= 1)
   {stop("For alternative 'two.sided': delta must be a positive value between 0 and p.hyp AND 0 and 1-p.hyp")} 

  }



# P.Ind gives TRUE or FALSE, depending on whether the interval includes the hypothetical proportion or not

 P.Ind.bin<-function(y,n,p.hyp,conf.level,alternative,method)
 {
  if(method=="Wald"){int=binWald(y=y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="Score"){int=binWilson(y=y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="AC"){int=binAC(y=y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="SOC"){int=binSOC(y=y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="CP"){int=binCP(y=y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="Blaker"){int=binBlaker(y=y,n=n, conf.level=conf.level)}

  return(int[1]>=p.hyp || int[2]<=p.hyp)
 }

# Probability of a certain event

 bin.prob <- function(y,n,p)
  {
   exp( lchoose(n,y) + y*log(p) + (n-y)*log(1-p) )
  }

yvec<-0:n

 if(alternative=="less")
  {
   p.tr=p.hyp-delta
   powvec<-numeric(length=length(yvec))   
   probvec<-numeric(length=length(yvec))

   for(i in 1:length(yvec))
    {powvec[i]<-P.Ind.bin(y=yvec[i], n=n, p.hyp=p.hyp, conf.level=conf.level, alternative=alternative, method=method)
     probvec[i]<-bin.prob(y=yvec[i], n=n, p=p.tr)
    }
  power=sum(powvec * probvec)

  }

 if(alternative=="greater")
  {
   p.tr=p.hyp+delta
   powvec<-numeric(length=length(yvec))   
   probvec<-numeric(length=length(yvec))

   for(i in 1:length(yvec))
    {powvec[i]<-P.Ind.bin(y=yvec[i], n=n, p.hyp=p.hyp, conf.level=conf.level, alternative=alternative, method=method)
     probvec[i]<-bin.prob(y=yvec[i], n=n, p=p.tr)
    }
  power=sum(powvec * probvec)
  }

 if(alternative=="two.sided")
  {
   p.trl <- p.hyp-delta
   p.tru <- p.hyp+delta
   powvec <- numeric(length=length(yvec))   
   probvecl <- numeric(length=length(yvec))
   probvecu <- numeric(length=length(yvec))

   for(i in 1:length(yvec))
    {powvec[i] <- P.Ind.bin(y=yvec[i], n=n, p.hyp=p.hyp, conf.level=conf.level, alternative=alternative, method=method)
     probvecl[i] <- bin.prob(y=yvec[i], n=n, p=p.trl)
     probvecu[i] <- bin.prob(y=yvec[i], n=n, p=p.tru)
    }

   powerl <- sum(powvec * probvecl)
   poweru <- sum(powvec * probvecu)
   power <- min(poweru,powerl)
  }

out<-list(power=power)

return(out)
}

