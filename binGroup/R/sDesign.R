"sDesign" <-
function(n,smax,delta,p.hyp,conf.level=0.95, power=0.8, alternative="two.sided", method="CP", biasrest=0.05)

{

 if(length(smax)!=1 || (smax<3 | abs(round(smax)-smax) > 1e-07))
  {stop("the maximal group size smax allowed in calculations must be a single integer greater than 0")}
 if(length(n)!=1 || (n<=1 | abs(round(n)-n) > 1e-07))
  {stop("the number of groups n must be specified as a single integer>1")}
 if(length(conf.level)!=1 || conf.level<0 || conf.level>1)
  {stop("conf.level must be a positive number between 0 and 1")}
 if(length(power)!=1 || power<0 || power>1)
  {stop("desired power must be a positive number between 0 and 1, f.e. 0.8 for rejecting H0 in 80% of the cases")}

  method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))

  alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))
 
 if(length(p.hyp)!=1 || p.hyp>1 || p.hyp<0) 
  {stop("threshold p.hyp must be specified as a single number between 0 and 1")}

 if( length(delta)!=1)
  {stop("delta must be specified as a single number")}
 if(alternative=="less")
  {
  if( p.hyp-delta < 0 || p.hyp-delta > 1 )
   {stop("alternative=less: specify delta as a number between 0 and the threshold p.hyp")}
  }

 if(alternative=="greater")
  {
  if( p.hyp+delta < 0 || p.hyp+delta > 1 )
   {stop("alternative=greater: specify delta as a number between the threshold p.hyp and 1")}
  }

 if(alternative=="two.sided")
  {
  if( p.hyp+delta < 0 || p.hyp+delta > 1 || p.hyp-delta < 0 || p.hyp-delta > 1)
   {stop("alternative=two.sided: specify delta as a number between the threshold p.hyp and 1")} 
  }

 if(length(biasrest)!=1 || biasrest>=1 || biasrest<0)
  {stop("the maximally allowed bias(p) specified in biasrest must be a single number between 0 and 1, usually should be close to 0")}


# # # Iteration until smax, until either the desired power is reached or biasrestriction is violated 


if(method=="SOC" && n<=3){stop("number of groups n<=3 might cause problems in computation of SOC interval")}


 sit<-2:smax
 powerit<-numeric(length=length(sit))
 biasit<-numeric(length=length(sit))

for (i in 1:length(sit))
  {

  temp=bgtPowerI(n=n, s=sit[i], delta=delta, p.hyp=p.hyp, conf.level=conf.level, alternative=alternative, method=method)
  powerit[i]<-temp$power
  biasit[i]<-temp$bias

  if(temp$bias <= biasrest && temp$power >= power)
    {
    out <- list(sout=sit[i], powerout=powerit[i], biasout=biasit[i], 
     power.reached=TRUE, bias.reached=FALSE, powerit=powerit, biasit=biasit, sit=sit, maxit=i,
     alternative=alternative, p.hyp=p.hyp, delta=delta, biasrest=biasrest, power=power)

     class(out)<-"sDesign"
     return(out)
    }

  if(temp$bias > biasrest)
    {
    out <- list(sout=sit[which.max(powerit[1:(i-1)])],
     powerout=powerit[which.max(powerit[1:(i-1)])],
     biasout=biasit[which.max(powerit[1:(i-1)])],
     power.reached=FALSE, bias.reached=TRUE,  
     powerit=powerit,biasit=biasit,sit=sit, maxit=i,
     alternative=alternative, p.hyp=p.hyp, delta=delta, biasrest=biasrest, power=power)

     class(out)<-"sDesign"
     return(out)
    }  
  }
## end of for statement

     out <- list(sout=sit[which.max(powerit)],
       powerout=powerit[which.max(powerit)],
       biasout=biasit[which.max(powerit)], 
       power.reached=FALSE, bias.reached=FALSE,
       powerit=powerit, biasit=biasit,sit=sit,maxit=length(sit),
       alternative=alternative, p.hyp=p.hyp, delta=delta, biasrest=biasrest, power=power )

     class(out)<-"sDesign"
     return(out)
}

