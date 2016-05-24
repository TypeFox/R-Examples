"binDesign" <-
function(nmax, delta, p.hyp, conf.level=0.95, power=0.8, method="CP", alternative="two.sided")
{

if( min(nmax) < 4 ) {stop(" nmax must be at least 4 ")} 
if(length(nmax)<1 || length(nmax)>2 ) {stop(" nmax must be either a single integer or a vector of two integers")} 

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



 if(length(nmax)==2){nit<-min(nmax):max(nmax)}
 if(length(nmax)==1){nit<-3:nmax}

 powerit=numeric(length=length(nit))
  for(n in 1:length(nit))
   {
    powerit[n] <- binPower(n=nit[n], delta=delta, p.hyp=p.hyp, conf.level=conf.level,
              method=method, alternative=alternative)$power
    if(powerit[n]>=power)     
     {out <- list(powerout=powerit[n],
               nout=nit[n],
               power.reached=TRUE,
               powerit=powerit,nit=nit,maxit=n,
               delta=delta, p.hyp=p.hyp, power=power, method=method, alternative=alternative)

     class(out) <- "binDesign"     
     return(out)
     }
  }

 npowmax <- nit[which.max(powerit)]
 powout <- powerit[which.max(powerit)] 
 out <- list(powerout=powout,
              nout=npowmax,
              power.reached=FALSE,
              powerit=powerit, nit=nit, maxit=length(nit),
              delta=delta, p.hyp=p.hyp, power=power, method=method, alternative=alternative)
 class(out) <- "binDesign"
 return(out)

}

