# We built a method for computing the density of Mixed Tempered Stable
setGeneric("dMixedTS",
            function(object,x=numeric(),setSup=NULL,setInf=NULL,N=2^10, ...)
              standardGeneric("dMixedTS")
)

setMethod("dMixedTS","param.MixedTS",
          function(object,x=numeric(),setSup=NULL,setInf=NULL,N=2^10,...){
            if(object@Mixing=="Gamma"){
              density<-dMixedTS.aux(x=x, mu0=object@mu0,
                       mu=object@mu, sig=object@sigma, a=object@a,
                       alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                       setSup=setSup, setInf=setInf, MixingDens=object@Mixing,N=N,
                       UseMGF=NULL, paramMixing=NULL, Parametrization=object@Parametrization)
              resDens<-new("MixedTS", dens=density, xMixedTS=x, quantile=FALSE, mu0=object@mu0,
                           mu=object@mu, sig=object@sigma, a=object@a,
                           alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                           Mixing=object@Mixing, Parametrization=object@Parametrization)
              
            }else{
              density<-dMixedTS.aux(x=x, mu0=object@mu0,
                                    mu=object@mu, sig=object@sigma, a=object@a,
                                    alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                                    setSup=setSup, setInf=setInf, MixingDens=object@Mixing,
                                    N=N,UseMGF=object@MixingLogMGF, paramMixing=object@paramMixing,
                                    Parametrization=object@Parametrization)
              resDens<-new("MixedTS", dens=density, xMixedTS=x, quantile=FALSE, mu0=object@mu0,
                           mu=object@mu, sig=object@sigma, a=object@a,
                           alpha=object@alpha, lambda_p=object@lambda_p,lambda_m=object@lambda_m,
                           Mixing=object@Mixing,MixingLogMGF=object@MixingLogMGF,
                           paramMixing=object@paramMixing,Parametrization=object@Parametrization)
            }            
            return(resDens)
          }
)

setGeneric("plot")


setMethod("plot",signature(x="MixedTS"),
          function(x,Nclass=20,type="l", main="Dens.MixedTS",ylab="Dens",...){
            obs<-NULL
            yvar<-NULL
            if(length(x@Data)==0){
              if(length(x@dens)>0){
                yvar<-x@dens
                obs<-x@xMixedTS
              }
              if(length(x@prob)>0 && x@quantile==FALSE){
                yvar<-x@prob
                obs<-x@xMixedTS
                if(main=="Dens.MixedTS"){
                  main<-"Prob. MixedTS"
                }
              }
              if(x@quantile==TRUE){
                yvar<-x@xMixedTS
                obs<-x@prob
                if(main=="Dens.MixedTS"){
                  main<-"Quantile MixedTS"
                }
              }
                .local<-function(obs,yvar,main,x,type,Nclass,...){
                  plot(x=obs,y=yvar,type=type, main=main,...)
                }
                
               }else{
                 #obs<-seq(min(x@Data),max(x@Data), length=Nclass)
                 obs<-NULL
                 yvar<-x@Data
                 if(main=="Dens.MixedTS"){
                   main<-"Rand. MixedTS"
                 }
                 .local<-function(obs,yvar,main,x,type,Nclass,...){
                       obs<-yvar
                        hist.Data<-hist(obs,nclass=Nclass,freq=FALSE,main=main, ...)
                       densdum<-dMixedTS(object=x,x=hist.Data$mids, N=2^12)
#                       plot(x=hist.Data$mids,y=hist.Data$dens,type="h")
                       lines(x=hist.Data$mids, y=densdum@dens, col="red")
                 }
            }
          .local(obs,yvar,main,x,type, Nclass,...)
          }
)
