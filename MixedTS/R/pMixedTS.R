# This function is for building the cdf of a Mixed TS. Internally the 
# function uses the pMixedTS.aux

setGeneric("pMixedTS",
           function(object,x,setSup=NULL,setInf=NULL,N=2^10, ...)
             standardGeneric("pMixedTS")
)

setMethod("pMixedTS","param.MixedTS",
          function(object,x,setSup=NULL,setInf=NULL,N=2^10,...){
            if(object@Mixing=="Gamma"){
              prob<-pMixedTS.aux(q=x, mu0=object@mu0,
                                    mu=object@mu, sig=object@sigma, a=object@a,
                                    alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                                    setSup=setSup, setInf=setInf, MixingDens=object@Mixing,N=N,
                                    UseMGF=NULL, paramMixing=NULL, 
                                    Parametrization=object@Parametrization)
              resDens<-new("MixedTS", prob=prob, xMixedTS=x, quantile=FALSE, mu0=object@mu0,
                           mu=object@mu, sig=object@sigma, a=object@a,
                           alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                           Mixing=object@Mixing,Parametrization=object@Parametrization)
              
            }else{
              prob<-pMixedTS.aux(q=x, mu0=object@mu0,
                                    mu=object@mu, sig=object@sigma, a=object@a,
                                    alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                                    setSup=setSup, setInf=setInf, MixingDens=object@Mixing,
                                    N=N,UseMGF=object@MixingLogMGF, 
                                 paramMixing=object@paramMixing, 
                                 Parametrization=object@Parametrization)
              resDens<-new("MixedTS", prob=prob, xMixedTS=x, quantile=FALSE, mu0=object@mu0,
                           mu=object@mu, sig=object@sigma, a=object@a,
                           alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                           Mixing=object@Mixing,MixingLogMGF=object@MixingLogMGF, 
                           paramMixing=object@paramMixing,
                           Parametrization=object@Parametrization)
            }            
            return(resDens)
          }
)
