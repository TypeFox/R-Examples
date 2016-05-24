# We built a method for computing the density of Mixed Tempered Stable
setGeneric("rMixedTS",
           function(object,x,setSup=10,setInf=-10,N=2^10, ...)
             standardGeneric("rMixedTS")
)

setMethod("rMixedTS","param.MixedTS",
          function(object,x,setSup=10,setInf=-10,N=2^10,...){
            if(object@Mixing=="Gamma"){
              ran<-rMixedTS.aux(n=x, mu0=object@mu0,
                                    mu=object@mu, sig=object@sigma, a=object@a,
                                    alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                                    setSup=setSup, setInf=setInf, MixingDens=object@Mixing,Nstep=N,
                                    UseMGF=NULL, paramMixing=NULL,
                                Parametrization=object@Parametrization)
              resDens<-new("MixedTS", Data=ran, xMixedTS=x, quantile=FALSE, mu0=object@mu0,
                           mu=object@mu, sig=object@sigma, a=object@a,
                           alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                           Mixing=object@Mixing,Parametrization=object@Parametrization)
              
            }else{
              ran<-rMixedTS.aux(n=x, mu0=object@mu0,
                                    mu=object@mu, sig=object@sigma, a=object@a,
                                    alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                                    setSup=setSup, setInf=setInf, MixingDens=object@Mixing,
                                    Nstep=N,UseMGF=object@MixingLogMGF, 
                                paramMixing=object@paramMixing,Parametrization=object@Parametrization)
              resDens<-new("MixedTS", Data=ran, xMixedTS=x, quantile=FALSE, mu0=object@mu0,
                           mu=object@mu, sig=object@sigma, a=object@a,
                           alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                           Mixing=object@Mixing,MixingLogMGF=object@MixingLogMGF,
                           paramMixing=object@paramMixing,Parametrization=object@Parametrization)
            }            
            return(resDens)
          }
)