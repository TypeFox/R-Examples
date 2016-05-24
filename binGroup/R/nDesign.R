"nDesign" <-
function(nmax,s,delta,p.hyp,conf.level=0.95, power=0.8, alternative="two.sided", method="CP", biasrest=0.05)
{

 if(length(nmax)<1 || length(nmax)>2 || ( min(nmax)<=3 | abs(round(nmax)-nmax) > 1e-07))
  {stop("the maximal number of groups n allowed in calculations must be one or two integer(s) greater than 1")}
 if(length(s)!=1 || (s<1 | abs(round(s)-s) > 1e-07))
  {stop("group size s must be specified as a single integer>0")}
 if(length(conf.level)!=1 || conf.level<0 || conf.level>1)
  {stop("conf.level must be a positive number between 0 and 1")}
 if(length(power)!=1 || power<0 || power>1)
  {stop(" desired power must be a positive number between 0 and 1, f.e. 0.8 for rejecting H0 in 80% of the cases")}

 method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))

 alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

 if(length(p.hyp)!=1 || p.hyp>1 || p.hyp<0)
  {stop("true proportion p.hyp must be specified as a single number between 0 and 1")}

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

# # # # # # 


if(length(nmax)==1)
{
 nit<-4:nmax
 powerit<-numeric(length=length(nit))
 biasit<-numeric(length=length(nit))

 for (i in 1:length(nit))
  {
  temp<-bgtPowerI(n=nit[i], s=s, delta=delta, p.hyp=p.hyp, conf.level=conf.level, alternative=alternative, method=method)
  powerit[i]<-temp$power
  biasit[i]<-temp$bias
  if(temp$bias <= biasrest && temp$power >= power)
    {   
        out<-list(nout=nit[i], powerout=powerit[i], biasout=temp$bias, power.reached=TRUE, bias.reached=FALSE, biasit=biasit,
        nit=nit, powerit=powerit, delta=delta, p.hyp=p.hyp, power=power, biasrest=biasrest, alternative=alternative, maxit=i ) 
        class(out)<-"nDesign"
        return(out)
    }
  }

# # # bias decreases monotone with increasing n: if nmax has bias> biasrest: all designs have bias> biasrest
lastn<-length(nit)
 if(biasit[lastn] > biasrest)
  { 
    
        out<-list(nout=nmax, powerout=powerit[lastn], biasout=biasit[lastn], power.reached=FALSE, bias.reached=TRUE,biasit=biasit,
        nit=nit, powerit=powerit,delta=delta,p.hyp=p.hyp, power=power, biasrest=biasrest, alternative=alternative, maxit=lastn) 
        class(out)<-"nDesign"
        return(out)  
  }

 npowmax=nit[which.max(powerit)]
 powout=powerit[which.max(powerit)] 
 biasout=biasit[which.max(powerit)]
  {
      out<-list(nout=npowmax, powerout=powout, biasout=biasout, power.reached=FALSE, bias.reached=FALSE,biasit=biasit,
        nit=nit, powerit=powerit,delta=delta,p.hyp=p.hyp, power=power, biasrest=biasrest, alternative=alternative,         maxit=length(nit)) 
        class(out)<-"nDesign"
        return(out)  

  }

}



if(length(nmax)==2){

 nfrom=min(nmax)
 nto=max(nmax)
 if(nfrom<4 && method=="SOC"){stop("The SOC interval can have bounds NaN for n < 4")}

 nit<-nfrom:nto
 powerit<-numeric(length=length(nit))
 biasit<-numeric(length=length(nit))

 for (i in 1:length(nit))
  {
  temp<-bgtPowerI(n=nit[i], s=s, delta=delta, p.hyp=p.hyp, conf.level=conf.level, alternative=alternative, method=method)
  powerit[i]<-temp$power
  biasit[i]<-temp$bias
  if(temp$bias <= biasrest && temp$power >= power)
    {   
        out<-list(nout=nit[i], powerout=powerit[i], biasout=temp$bias, power.reached=TRUE, bias.reached=FALSE,biasit=biasit,
        nit=nit, powerit=powerit,delta=delta,p.hyp=p.hyp, power=power, biasrest=biasrest, alternative=alternative, maxit=i )
        class(out)<-"nDesign"
        return(out)
    }
  }

# # # bias decreases monotone with increasing n: if nmax has bias> biasrest: all designs have bias> biasrest
lastn<-length(nit)
 if(biasit[lastn] > biasrest)
  {
        out<-list(nout=max(nmax), powerout=powerit[lastn], biasout=biasit[lastn], power.reached=FALSE,biasit=biasit,                 bias.reached=TRUE,nit=nit, powerit=powerit,delta=delta,p.hyp=p.hyp, power=power, biasrest=biasrest,                 alternative=alternative, maxit=lastn) 
        class(out)<-"nDesign"
        return(out)   
  }

 npowmax=nit[which.max(powerit)]
 powout=powerit[which.max(powerit)] 
 biasout=biasit[which.max(powerit)]
  {    
        out<-list(nout=npowmax, powerout=powout, biasout=biasout, power.reached=FALSE, bias.reached=FALSE,biasit=biasit,
        nit=nit, powerit=powerit,delta=delta,p.hyp=p.hyp, power=power, biasrest=biasrest, alternative=alternative, maxit=length(nit)) 
        class(out)<-"nDesign"
        return(out)

  }


}




}

