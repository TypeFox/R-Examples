# Main Classes in the MixedTS package.
# Four main classes
# I) parameters Mixed Tempered Stable Classes that allows the user to specify the set
# of parameters.
setClass("param.MixedTS",representation(mu0="numeric",
                                 mu="numeric",
                                 sigma="numeric",
                                 a="vector",
                                 alpha="numeric",
                                 lambda_p="numeric",
                                 lambda_m="numeric",
                                 Mixing="character",
                                 paramMixing="list",
                                 MixingLogMGF="OptionalFunction",
                                 Parametrization = "character"
                                 ),
                          prototype(mu0=numeric(),
                                    mu=numeric(),
                                    sigma=numeric(),
                                    alpha=numeric(),
                                    lambda_p=numeric(),
                                    lambda_m=numeric(),
                                    Mixing=character(),
                                    paramMixing=list(),
                                    MixingLogMGF=NULL,
                                    Parametrization = character())
          )
# 
setMethod("initialize", "param.MixedTS",
          function(.Object,
                   mu0=numeric(),
                   mu=numeric(),
                   sigma=numeric(),
                   a=numeric(),
                   alpha=numeric(),
                   lambda_p=numeric(),
                   lambda_m=numeric(),
                   Mixing=character(),
                   paramMixing=list(),
                   MixingLogMGF=function(){},
                   Parametrization=character()
                   ){
            .Object@mu0 <- mu0
            .Object@mu <- mu
            .Object@sigma <- sigma
            .Object@a <- a
            .Object@alpha <- alpha
            .Object@lambda_p <- lambda_p
            .Object@lambda_m <- lambda_m
            .Object@Mixing <- Mixing
            .Object@paramMixing <- paramMixing
            .Object@MixingLogMGF<-MixingLogMGF
            .Object@Parametrization <-Parametrization
            .Object
          }
          )
# This function is used internally  by setMixedTS and returns a object of class 
# param.MixedTS.
setMixedTS.param <-function(mu0=numeric(),
                            mu=numeric(),
                            sigma=numeric(),
                            a,
                            alpha=numeric(),
                            lambda_p=numeric(),
                            lambda_m=numeric(),
                            param=numeric(),
                            Mixing="Gamma",
                            paramMixing=list(),
                            Parametrization="A"
                            ){
  if(length(param)==0){
    if(sigma<0){
      warning("sigma must be positive")
      return(NULL)
    }
    if(Mixing=="Gamma"){
    if(a<0){
      warning("a must be positive")
      return(NULL)      
    }
    MixingLogMGF=NULL
    }else{
      if(Mixing=="User"){
        if(!is.character(a)){
          warning("a must be a string that is the formula of the Log MGF of Mixing Random Variable")
          return(NULL)
          }else{
          if(!is.list(paramMixing)){
            warning("The mixing parameters must be stored in a list")
            return(NULL)
          }
          charactFunc<-function(u,paramMixing,a){
            u<-u
            N.Inp<-names(paramMixing)
            for(i in 1: length(N.Inp)){
              assign(N.Inp[i],as.complex(paramMixing[i]))
            }
            eval(parse(text=a))
          }
          MixingLogMGF<-charactFunc
        }
      }  
    }
    if((alpha<0)&&(alpha>2)){
      warning("alpha is outside of its domain")
      return(NULL)
    }
    if(lambda_p<0){
      warning("lambda_p must be positive")
      return(NULL)
    }
    if(lambda_m<0){
      warning("lambda_m must be positive")
      return(NULL)
      
    }
    if(Parametrization!="A"){
      if(Parametrization!="B"){
        warning("Choose A or B for Paramtrization arg. See Help for the meaning")
        return(NULL)
      }
    }
    res<-new("param.MixedTS",
             mu0=mu0,
             mu=mu,
             sigma=sigma,
             a=a,
             alpha=alpha,
             lambda_p=lambda_p,
             lambda_m=lambda_m,
             Mixing=Mixing,
             paramMixing=paramMixing,
             MixingLogMGF=MixingLogMGF,
             Parametrization=Parametrization)
#     res@Mixing<-Mixing
#     res@paramMixing<-paramMixing
#     res@MixingLogMGF<-MixingLogMGF
    return(res)
  }else{
    if(param[3]<0){
      warning("sigma must be positive")
      return(NULL)
    }
    if(param[4]<0){
      warning("a must be positive")
      return(NULL)
      
    }
    if((param[5]<0)&&(param[5]>2)){
      warning("alpha is outside of its domain")
      return(NULL)
    }
    if(param[6]<0){
      warning("lambda_p must be positive")
      return(NULL)
    }
    if(param[7]<0){
      warning("lambda_m must be positive")
      return(NULL)
      
    }
    
    res<-new("param.MixedTS",
             mu0=param[1],
             mu=param[2],
             sigma=param[3],
             a=param[4],
             alpha=param[5],
             lambda_p=param[6],
             lambda_m=param[7],
             Mixing=Mixing,
             paramMixing=paramMixing,
             MixingLogMGF=MixingLogMGF,
             Parametrization=Parametrization)
    
  }
}
# Main class of the package MixedTS. This Class contains all info relating to
# the distribution and eventaully the data.
MixedTSClass <- setClass("MixedTS", 
                       representation(Data="numeric",
                                      dens="numeric",
                                      prob="numeric",
                                      xMixedTS="numeric",
                                      quantile="logical"),
                       prototype(Data=numeric(),
                                 dens=numeric(),
                                 prob=numeric(),
                                 xMixedTS=numeric(),
                                 quantile=FALSE),
                       contains="param.MixedTS"
)
# 
 setMethod("initialize", "MixedTS",
             function(.Object,
                      Data=numeric(),
                      dens=numeric(),
                      prob=numeric(),
                      xMixedTS=numeric(),
                      quantile=logical(),
                      mu0=numeric(),
                      mu=numeric(),
                      sigma=numeric(),
                      a=numeric(),
                      alpha=numeric(),
                      lambda_p=numeric(),
                      lambda_m=numeric(),
                      Mixing=character(),
                      paramMixing=list(),
                      MixingLogMGF=function(){},
                      Parametrization=character()){
               .Object@Data=Data
               .Object@Mixing=Mixing
               .Object@dens=dens
               .Object@prob=prob
               .Object@xMixedTS=xMixedTS
               .Object@mu0 <- mu0
               .Object@mu <- mu
               .Object@sigma <- sigma
               .Object@a <- a
               .Object@alpha <- alpha
               .Object@lambda_p <- lambda_p
               .Object@lambda_m <- lambda_m
               .Object@Mixing <- Mixing
               .Object@paramMixing <- paramMixing
               .Object@MixingLogMGF<-MixingLogMGF
               .Object@quantile<-quantile
               .Object@Parametrization<-Parametrization
               .Object
              },
                      
)
# 
# # setMixedTS<-function(param=numeric(),
# #                      mu0=numeric(),
# #                      mu=numeric(),
# #                      a=numeric(),
# #                      alpha=numeric(),
# #                      lambda_p=numeric(),
# #                      lambda_m=numeric(),
# #                      Mixing="Gamma",
# #                      Data=numeric(),
# #                      xMin=numeric(),
# #                      xMax=numeric(),
# #                      nStep=100,
# #                      qMin=0.01,
# #                      qMax=0.99){
# # 
# #                      ParamMTS <-setMixedTS.param(mu0=mu0,
# #                                             mu=mu,
# #                                             sigma=sigma,
# #                                             a=a,
# #                                             alpha=alpha,
# #                                             lambda_p=lambda_p,
# #                                             lambda_m=lambda_m,
# #                                             param=param)
# #                      if(length(Data)==0){
# #                        if((length(xMin)>0)&&(length(xMin)>0)){
# #                          xMTS<-seq(xMin,xMax,length=nStep)
# #                        }else{
# #                          xMTS<-qMixedTS.aux(p=seq(qMin,qMax,length=nStep),
# #                                   mu0=ParamMTS@mu0,
# #                                   mu=ParamMTS@mu,
# #                                   sig=ParamMTS@sigma,
# #                                   a=ParamMTS@a,
# #                                   alpha=ParamMTS@alpha,
# #                                   lambda_p=ParamMTS@lambda_p,
# #                                   lambda_m=ParamMTS@lambda_m)
# #                        }
# #                      }else{
# #                         xMTS<-seq(min(Data),max(Data),length=nStep)
# #                       }
# #                      
# #                      densMTS<-dMixedTS.aux(x=xMTS,
# #                                        mu0=ParamMTS@mu0,
# #                                        mu=ParamMTS@mu,
# #                                        sig=ParamMTS@sigma,
# #                                        a=ParamMTS@a,
# #                                        alpha=ParamMTS@alpha,
# #                                        lambda_p=ParamMTS@lambda_p,
# #                                        lambda_m=ParamMTS@lambda_m)
# #                      ProbMTS<-pMixedTS.aux(q=xMTS,
# #                                        mu0=ParamMTS@mu0,
# #                                        mu=ParamMTS@mu,
# #                                        sig=ParamMTS@sigma,
# #                                        a=ParamMTS@a,
# #                                        alpha=ParamMTS@alpha,
# #                                        lambda_p=ParamMTS@lambda_p,
# #                                        lambda_m=ParamMTS@lambda_m)
# #                      ResMTS<-new("MixedTS",
# #                                 param=ParamMTS,
# #                                 Data=Data,
# #                                 Mixing=Mixing,
# #                                 dens=densMTS,
# #                                 prob=ProbMTS,
# #                                 xMixedTS=xMTS)
# #                      return(ResMTS)
# #   
# # }
# #
# 
# setMixedTS<-function(param=numeric(),
#                      mu0=numeric(),
#                      mu=numeric(),
#                      a=numeric(),
#                      alpha=numeric(),
#                      lambda_p=numeric(),
#                      lambda_m=numeric(),
#                      Mixing="Gamma",
#                      Data=numeric(),
#                      xMin=numeric(),
#                      xMax=numeric(),
#                      nStep=100,
#                      qMin=0.01,
#                      qMax=0.99){
#   
#   ParamMTS <-setMixedTS.param(mu0=mu0,
#                               mu=mu,
#                               sigma=sigma,
#                               a=a,
#                               alpha=alpha,
#                               lambda_p=lambda_p,
#                               lambda_m=lambda_m,
#                               param=param)
#   if(length(Data)==0){
#     if((length(xMin)>0)&&(length(xMin)>0)){
#       xMTS<-seq(xMin,xMax,length=nStep)
#     }else{
#       xMTS<-qMixedTS.aux(p=seq(qMin,qMax,length=nStep),
#                      mu0=ParamMTS@mu0,
#                      mu=ParamMTS@mu,
#                      sig=ParamMTS@sigma,
#                      a=ParamMTS@a,
#                      alpha=ParamMTS@alpha,
#                      lambda_p=ParamMTS@lambda_p,
#                      lambda_m=ParamMTS@lambda_m)
#     }
#   }else{
#     xMTS<-seq(min(Data),max(Data),length=nStep)
#   }
#   
#   densMTS<-dMixedTS.aux(x=xMTS,
#                     mu0=ParamMTS@mu0,
#                     mu=ParamMTS@mu,
#                     sig=ParamMTS@sigma,
#                     a=ParamMTS@a,
#                     alpha=ParamMTS@alpha,
#                     lambda_p=ParamMTS@lambda_p,
#                     lambda_m=ParamMTS@lambda_m)
#   ProbMTS<-pMixedTS.aux(q=xMTS,
#                     mu0=ParamMTS@mu0,
#                     mu=ParamMTS@mu,
#                     sig=ParamMTS@sigma,
#                     a=ParamMTS@a,
#                     alpha=ParamMTS@alpha,
#                     lambda_p=ParamMTS@lambda_p,
#                     lambda_m=ParamMTS@lambda_m)
#   ResMTS<-new("MixedTS",
#               param=ParamMTS,
#               Data=Data,
#               Mixing=Mixing,
#               dens=densMTS,
#               prob=ProbMTS,
#               xMixedTS=xMTS)
#   return(ResMTS)
#   
# }
# 


  MixedTSqmle <- setClass("MixedTS.qmle",
                          representation(time="numeric",
                                         coef="numeric",
                                         vcov="matrix",
                                         min="numeric",
                                         details="list",
                                         nobs="integer",
                                         method="character"),
                          prototype(time=numeric(),
                                    coef=numeric(),
                                    vcov=matrix(),
                                    min=numeric(),
                                    details=list(),
                                    nobs=integer(),
                                    method=character()
                                    ),
                          contains=c("MixedTS")
 )
# # 
# 

setMethod("initialize", "MixedTS.qmle",
          function(.Object,
                   time=numeric(),
                   Data=numeric(),
                   dens=numeric(),
                   prob=numeric(),
                   xMixedTS=numeric(),
                   quantile=logical(),
                   mu0=numeric(),
                   mu=numeric(),
                   sigma=numeric(),
                   a=vector(),
                   alpha=numeric(),
                   lambda_p=numeric(),
                   lambda_m=numeric(),
                   Mixing=character(),
                   paramMixing=list(),
                   MixingLogMGF=function(){},
                   Parametrization=character(),
                   coef=numeric(),
                   vcov=matrix(),
                   min=numeric(),
                   details=list(),
                   nobs=integer(),
                   method=character()
          ){
            .Object@time=time
            .Object@Data=Data
            .Object@Mixing=Mixing
            .Object@dens=dens
            .Object@prob=prob
            .Object@xMixedTS=xMixedTS
            .Object@mu0 <- mu0
            .Object@mu <- mu
            .Object@sigma <- sigma
            .Object@a <- a
            .Object@alpha <- alpha
            .Object@lambda_p <- lambda_p
            .Object@lambda_m <- lambda_m
            .Object@Mixing <- Mixing
            .Object@paramMixing <- paramMixing
            .Object@MixingLogMGF<-MixingLogMGF
            .Object@quantile<-quantile
            .Object@coef<-coef
            .Object@vcov<-vcov
            .Object@min<-min
            .Object@details<-details
            .Object@nobs<-nobs
            .Object@method<-method
            .Object
          }
)
