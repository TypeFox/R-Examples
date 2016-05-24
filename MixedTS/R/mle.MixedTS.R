# # setGeneric("confint",
# #             function(object,...){
# #               standardGeneric("confint")
# #             }
# #            )
# # 
# # setMethod("confint","MixedTS.qmle",
# #            function(object,...){
# #              return(a="The function confint is not available for a class of object MixedTS.qmle")
# #            }
# #           )
# # 
# # setGeneric("profile",
# #            function(object,...){
# #              standardGeneric("profile")
# #            }
# # )
# # 
# # setMethod("profile","MixedTS.qmle",
# #           function(object, ...){
# #             return(a="The profile confint is not available for a class of object MixedTS.qmle")
# #           }
# # )
# 
# # We build a method for the Maximum Likelihood Estimation of The Mixed Tempered Stable.
# # This function use internally the function mle.MixedTS.aux in which if we fix the alpha=2,
# # we get the estimation  of the normal variance mean mixture.
# # The function is very flexible since allow us to estimate MixedTS with Gamma mixing 
# # density as well as the mixing density specified by the user.
# 
# setGeneric("mle.MixedTS",
#            function(object,Data=NULL,method="L-BFGS-B", fixed.param=NULL,
#                     lower.param=NULL,upper.param=NULL,setSup=NULL,setInf=NULL,N=2^10, ...)
#              standardGeneric("mle.MixedTS")
# )
# 
# setMethod("mle.MixedTS","param.MixedTS",
#           function(object,Data=NULL,method="L-BFGS-B", fixed.param=NULL,
#                    lower.param=NULL,upper.param=NULL,setSup=NULL,setInf=NULL,N=2^10,...){
#             call<-match.call()
#             
#             if(is(object,"MixedTS")){
#               if(length(object@Data)>0){
#                 if(is.null(Data)){
#                   Data<-object@Data
#                 }
#               }
#             }
#           
#           if(object@Mixing=="Gamma"){
#             param0<-list(mu0=object@mu0,mu=object@mu,sigma=object@sigma,a=object@a,
#                          alpha=object@alpha,lambda_p=object@lambda_p,lambda_m=object@lambda_m)
#             res<-mle.MixedTS.aux(data=Data, param0=param0, method, 
#                                  fixed.param, lower.param, upper.param,
#                                  MixingDens=object@Mixing, N, MixingLogMGF=NULL,
#                                  paramMixing=NULL,MGFdef=NULL,
#                                  setSup,setInf
#               )
# #            return(res)
#           }else{
#             param0<-list(mu0=object@mu0,mu=object@mu,sigma=object@sigma,
#                          alpha=object@alpha,lambda_p=object@lambda_p,lambda_m=object@lambda_m)
#             res<-mle.MixedTS.aux(data=Data, param0=param0, method, 
#                                  fixed.param, lower.param, upper.param,
#                                  MixingDens=object@Mixing, N, MixingLogMGF=object@MixingLogMGF,
#                                  paramMixing=object@paramMixing,MGFdef=object@a,
#                                  setSup,setInf
#             )
# #             dumm<-names(res$time)
# #             time<-as.numeric(res$time)
# #             names(time)<-dumm
# #             result<-new("MixedTS.qmle",time=time,Data=Data,
# #                     mu0=object@mu0, mu=object@mu,sigma=object@sigma,a=object@a,
# #                     alpha=object@alpha,lambda_p=object@lambda_p,lambda_m=object@lambda_m,
# #                     Mixing=object@Mixing, MixingLogMGF=object@MixingLogMGF,
# #                     paramMixing=object@paramMixing, 
# #                     coef=res$firs.prob$par,
# #                     fullcoef=res$firs.prob$par,
# #                     vcov=res$covErr,
# #                     min=res$firs.prob$value,
# #                     details=list(),
# #                     minuslogl=function(){},
# #                     nobs=length(Data),
# #                     method=method
# #                 )
#             
#           }
#           
#           
#           dumm<-names(res$time)
#           time<-as.numeric(res$time)
#           names(time)<-dumm
#           result<-new("MixedTS.qmle",time=time,Data=Data,
#                       mu0=object@mu0, mu=object@mu,sigma=object@sigma,a=object@a,
#                       alpha=object@alpha,lambda_p=object@lambda_p,lambda_m=object@lambda_m,
#                       Mixing=object@Mixing, MixingLogMGF=object@MixingLogMGF,
#                       paramMixing=object@paramMixing, 
#                       coef=res$firs.prob$par,
#                       call=call,
#                       fullcoef=res$firs.prob$par,
#                       vcov=res$covErr,
#                       min=res$firs.prob$value,
#                       details=res$firs.prob,
#                       minuslogl=minusloglik.MixedTS,
#                       nobs=length(Data),
#                       method=method
#           )
#           
#           return(result)
#   
#           
#   }
#   )

setGeneric("summary")
setMethod("summary","MixedTS.qmle",
          function(object,...){
          cmat <- cbind(Estimate = object@coef, `Std. Error` = sqrt(diag(object@vcov)))
          m2logL <- 2 * object@min
          new("summary.mle", call=call(quote("Mixed Tempered Stable")), coef = cmat, m2logL = m2logL)
}

          
          )
#             
setGeneric("coef")
setMethod("coef","MixedTS.qmle",
          function(object,...){
            object@coef 
          }
)

setGeneric("vcov")
setMethod("vcov","MixedTS.qmle",
          function(object,...){
            object@vcov 
          }
)
setGeneric("logLik")
setMethod("logLik","MixedTS.qmle",
          function(object,...){
            if (!missing(...)) 
              warning("extra arguments discarded")
            val <- -object@min
            if ("nobs" %in% slotNames(object) && !is.na(no <- object@nobs)) 
              attr(val, "nobs") <- no
            attr(val, "df") <- length(object@coef)
            class(val) <- "logLik"
            val
          }
)


setGeneric("BIC")
setMethod("BIC","MixedTS.qmle",
          function(object,...){
            logL <- logLik(object)
            npar<-length(object@coef)
            nobs <- object@nobs
            res<- -2*as.numeric(logL) + npar*log(nobs)
            names(res)<-"BIC"
            return(res)
          }
)

setGeneric("AIC")
setMethod("AIC","MixedTS.qmle",
          function(object,...){
            logL <- logLik(object)
            npar <- length(object@coef)
            nobs <- object@nobs
            res<-as.numeric(-2*logL + (2*nobs*npar)/(nobs-npar-1))
            names(res)<-"AIC"
            return(res)
          }
)



