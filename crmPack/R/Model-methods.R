#####################################################################################
## Author: Daniel Sabanes Bove[sabanesd *a*t* roche *.* com],
##         Wai Yin Yeung [w *.* yeung *a*t* lancaster *.* ac *.* uk]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Model-methods.R] by DSB Mon 11/05/2015 17:45>
##
## Description:
## Encapsulate the model input in a formal class.
##
## History:
## 31/03/2014   file creation
## 09/07/2015   Adding more methods for Pseudo Models
######################################################################################

##' @include Model-class.R
##' @include Samples-class.R
{}

## --------------------------------------------------
## Compute the doses for a given probability, given model and samples
## --------------------------------------------------

##' Compute the doses for a given probability, given model and samples
##'
##' @param prob the probability
##' @param model the \code{\linkS4class{Model}}
##' @param samples the \code{\linkS4class{Samples}}
##' @param \dots unused
##' 
##' @export
##' @keywords methods
setGeneric("dose",
           def=
           function(prob, model, samples, ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("dose")
           },
           valueClass="numeric")

##' @describeIn dose
##' @example examples/Model-method-dose.R
setMethod("dose",
          signature=
          signature(prob="numeric",
                    model="Model",
                    samples="Samples"),
          def=
          function(prob, model, samples, ...){
              ## extract the dose function from the model
              doseFun <- slot(model, "dose")
              ## which arguments, besides the prob, does it need?
              argNames <- setdiff(names(formals(doseFun)),
                                  "prob")
              ## now call the function with prob and with
              ## the arguments taken from the samples
              ret <- do.call(doseFun,
                             c(list(prob=prob),
                               samples@data[argNames]))
              ## return the resulting vector
              return(ret)
          })

## =====================================================================================
## -------------------------------------------------------------------------------------
## Compute the doses for a given probability, given Pseudo DLE model and given samples
## -----------------------------------------------------------------------------------
##' @describeIn dose Compute the doses for a given probability, given 
##' Pseudo DLE model with samples 
##' @example examples/Model-method-dose-modelTox.R
setMethod("dose",
          signature=
            signature(prob="numeric",
                      model="ModelTox",
                      samples="Samples"),
          def=
            function(prob, model, samples, ...){
              ## extract the dose function from the model
              doseFun <- slot(model, "dose")
              ## which arguments, besides the prob, does it need?
              argNames <- setdiff(names(formals(doseFun)),
                                  "prob")
              ## now call the function with prob and with
              ## the arguments taken from the samples
              ret <- do.call(doseFun,
                             c(list(prob=prob),
                               samples@data[argNames]))
              ## return the resulting vector
              return(ret)
            })

## =========================================================================
## ----------------------------------------------------------------------------
## Compute the dose for a given Pseudo DLE model and a given probability
## -----------------------------------------------------------------------
##' @describeIn dose Compute the dose for a given probability and a given 
##' Pseudo DLE model without samples
##' @example examples/Model-method-doseNoSamples.R
setMethod("dose",
          signature=
            signature(prob="numeric",
                      model="ModelTox",
                      samples="missing"),
          def=
            function(prob, model, ...){
              ## extract the dose function from the model
              doseFun <- slot(model, "dose")
              ## which arguments, besides the prob, does it need?
              argNames <- setdiff(names(formals(doseFun)),
                                  "prob")
              values<-c()
              for (parName in argNames){
                values<-c(values, slot(model,parName))}
              ## now call the function with prob
              ret <- do.call(doseFun,
                             c(list(prob=prob), values))
              ## return the resulting vector
              return(ret)
            })


## =======================================================================================

## --------------------------------------------------
## Compute the probability for a given dose, given model and samples
## --------------------------------------------------

##' Compute the probability for a given dose, given model and samples
##'
##' @param dose the dose
##' @param model the \code{\linkS4class{Model}} object
##' @param samples the \code{\linkS4class{Samples}}
##' @param \dots unused
##' @return the vector (for \code{\linkS4class{Model}} objects) of probability
##' samples.
##'
##' @export
##' @keywords methods
setGeneric("prob",
           def=
           function(dose, model, samples, ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("prob")
           })
## todo: simplify this to always vectorize internally over points -> always
## take/return matrix

##' @describeIn prob
##' @example examples/Model-method-prob.R
setMethod("prob",
          signature=
          signature(dose="numeric",
                    model="Model",
                    samples="Samples"),
          def=
          function(dose, model, samples, ...){
              ## extract the prob function from the model
              probFun <- slot(model, "prob")
              ## which arguments, besides the dose, does it need?
              argNames <- setdiff(names(formals(probFun)),
                                  "dose")
              ## now call the function with dose and with
              ## the arguments taken from the samples
              ret <- do.call(probFun,
                             c(list(dose=dose),
                               samples@data[argNames]))
              ## return the resulting vector
              return(ret)
          })

## =============================================================================================
## Compute the probability for a given dose, given Pseudo DLE model and samples
## --------------------------------------------------

##' @describeIn prob Compute the probability for a given dose, 
##' given Pseudo DLE model and samples
##' @example examples/Model-method-prob-modelTox.R
setMethod("prob",
          signature=
            signature(dose="numeric",
                      model="ModelTox",
                      samples="Samples"),
          def=
            function(dose, model, samples, ...){
              ## extract the prob function from the model
              probFun <- slot(model, "prob")
              ## which arguments, besides the dose, does it need?
              argNames <- setdiff(names(formals(probFun)),
                                  "dose")
              ## now call the function with dose and with
              ## the arguments taken from the samples
              ret <- do.call(probFun,
                             c(list(dose=dose),
                               samples@data[argNames]))
              ## return the resulting vector
              return(ret)
            })




## ==============================================================================
## ## Compute the probability for a given dose, given Pseudo DLE model
## --------------------------------------------------

##' @describeIn prob Compute the probability for a given dose, given Pseudo DLE model without samples
##' @example examples/Model-method-probNoSamples.R
setMethod("prob",
          signature=
            signature(dose="numeric",
                      model="ModelTox",
                      samples="missing"),
          def=
            function(dose,model,...){
              ## extract the prob function from the model
              probFun <- slot(model, "prob")
              ## which arguments, besides the dose, does it need?
              argNames <- setdiff(names(formals(probFun)),
                                  "dose")
              ## now call the function with dose
              values<-c()
              for (parName in argNames){
                values<-c(values, slot(model,parName))}
              ret <- do.call(probFun,
                             c(list(dose=dose), values))
              ## return the resulting vector
              return(ret)
            })

## =============================================================================

## --------------------------------------------------
## Compute the biomarker level for a given dose, given model and samples
## --------------------------------------------------

##' Compute the biomarker level for a given dose, given model and samples
##'
##' @param dose the dose
##' @param model the \code{\linkS4class{DualEndpoint}} object
##' @param samples the \code{\linkS4class{Samples}} object
##' @param \dots unused
##'  
##' @export
##' @keywords methods
setGeneric("biomLevel",
           def=
           function(dose, model, samples, ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("biomLevel")
           },
           valueClass="numeric")

##' @param xLevel the grid index of \code{dose}
##' @describeIn biomLevel Here it is very easy, we just return the corresponding
##' column (index \code{xLevel}) of the biomarker samples matrix, since we save
##' that in the samples
##' @example examples/Model-method-biomLevel.R
setMethod("biomLevel",
          signature=
          signature(dose="numeric",
                    model="DualEndpoint",
                    samples="Samples"),
          def=
          function(dose, model, samples, xLevel, ...){

              return(samples@data$betaW[, xLevel])

          })

## =============================================================================================
## ---------------------------------------------------------------------------------------------
## Compute the Expected Efficacy based on a given dose, a given pseduo Efficacy log-log model and a given 
## efficacy sample
## -----------------------------------------------------------------------------------------------
##' Compute the expected efficacy based on a given dose, a given pseudo Efficacy log-log model and a given 
##' efficacy sample
##' 
##' @param dose the dose
##' @param model the \code{\linkS4class{Effloglog}} class object
##' @param samples the \code{\linkS4class{Samples}} class object 
##' (can also be missing)
##' @param \dots unused
##' 
##' @example examples/Model-method-ExpEff.R
##' @export
##' @keywords methods
setGeneric("ExpEff",
           def=
             function(dose,model,samples,...){
               standardGeneric("ExpEff")
             },
           valueClass="numeric")

##' @describeIn ExpEff Method for the Effloglog class
setMethod("ExpEff",
          signature=
            signature(dose="numeric",
                      model="Effloglog",
                      samples="Samples"),
          def=
            function(dose, model, samples, ...){
              ## extract the ExpEff function from the model
              EffFun <- slot(model, "ExpEff")
              ## which arguments, besides the dose, does it need?
              argNames <- setdiff(names(formals(EffFun)),
                                  "dose")
              ## now call the function with dose and with
              ## the arguments taken from the samples
              ret <- do.call(EffFun,
                             c(list(dose=dose),
                               samples@data[argNames]))
              ## return the resulting vector
              return(ret)
            })
##======================================================================================

##' @describeIn ExpEff Compute the Expected Efficacy based a given dose and a given Pseudo Efficacy log log model without
##' samples
##' @example examples/Model-method-ExpEffNoSamples.R
setMethod("ExpEff",
          signature=
            signature(dose="numeric",
                      model="Effloglog",
                      samples="missing"),
          def=
            function(dose, model, ...){
              ## extract the ExpEff function from the model
              EffFun <- slot(model, "ExpEff")
              ## which arguments, besides the dose, does it need?
              argNames <- setdiff(names(formals(EffFun)),
                                  "dose")
              ## now call the function with dose
              values<-c()
              for (parName in argNames){
                values <- c(values, slot(model,parName))}
              
              ret <- do.call(EffFun,
                             c(list(dose=dose),values))
              ## return the resulting vector
              return(ret)
            })

##' @describeIn ExpEff Compute the Expected Efficacy based a given dose, Efficacy 
##' Flexible model with samples
##' @example examples/Model-method-ExpEffFlexi.R
setMethod("ExpEff",
          signature=
            signature(dose="numeric",
                      model="EffFlexi",
                      samples="Samples"),
          def=
            function(dose, model, samples, ...){
              ## extract the ExpEff function from the model
              EffFun <- slot(model, "ExpEff")
              ret <- EffFun(dose,data=model@data,Effsamples=samples)
              ## return the resulting vector
              return(ret)
            })

## ---------------------------------------------------------------------------------
## Compute gain value using a Pseudo DLE and a pseduo Efficacy log-log model
## -------------------------------------------------------------------------------

##' Compute the gain value with a given dose level, given a pseudo DLE model, a DLE sample, 
##' a pseudo Efficacy log-log model and a Efficacy sample
##' 
##' @param dose the dose
##' @param DLEmodel the \code{\linkS4class{ModelTox}} object
##' @param DLEsamples the \code{\linkS4class{Samples}} object (can also be missing)
##' @param Effmodel the \code{\linkS4class{Effloglog}} or the \code{\linkS4class{EffFlexi}} object
##' @param Effsamples the \code{\linkS4class{Samples}} object (can also be missing)
##' @param \dots unused
##' 
##' @export
##' @keywords methods
setGeneric("gain",
           def=
             function(dose,DLEmodel,DLEsamples,Effmodel,Effsamples,...){
               standardGeneric("gain")
             },
           valueClass="numeric")

##' @describeIn gain
##' @example examples/Model-method-gain.R
setMethod("gain",
          signature=
            signature(dose="numeric",
                      DLEmodel="ModelTox",
                      DLEsamples="Samples",
                      Effmodel="Effloglog",
                      Effsamples="Samples"),
          def=
            function(dose,DLEmodel,DLEsamples, Effmodel,Effsamples,...){
              
              
              ## extract the prob function from the model
              probFun <- slot(DLEmodel, "prob")
              ## which arguments, besides the dose, does it need?
              DLEargNames <- setdiff(names(formals(probFun)),
                                     "dose")
              ## now call the function with dose and with
              ## the arguments taken from the samples
              DLEret <- do.call(probFun,
                                c(list(dose=dose),
                                  DLEsamples@data[DLEargNames]))
              
              ## extract the ExpEff function from the model
              EffFun <- slot(Effmodel, "ExpEff")
              ## which arguments, besides the dose, does it need?
              EffargNames <- setdiff(names(formals(EffFun)),
                                     "dose")
              ## now call the function with dose and with
              ## the arguments taken from the samples
              Effret <- do.call(EffFun,
                                c(list(dose=dose),
                                  Effsamples@data[EffargNames]))
              
              ## return the resulting vector
              Gainret <- Effret/(1+(DLEret/(1-DLEret)))
              return(Gainret)
            })

## ===================================================================


##' @describeIn gain Compute the gain given a dose level, a pseduo DLE model, a DLE sample, 
##' the pseudo EffFlexi model and an Efficacy sample
##' @example examples/Model-method-gainFlexi.R
setMethod("gain",
          signature=
            signature(dose="numeric",
                      DLEmodel="ModelTox",
                      DLEsamples="Samples",
                      Effmodel="EffFlexi",
                      Effsamples="Samples"),
          def=
            function(dose,DLEmodel,DLEsamples, Effmodel,Effsamples,...){
              
              
              
              ## extract the prob function from the model
              probFun <- slot(DLEmodel, "prob")
              ## which arguments, besides the dose, does it need?
              DLEargNames <- setdiff(names(formals(probFun)),
                                     "dose")
              ## now call the function with dose and with
              ## the arguments taken from the samples
              DLEret <- do.call(probFun,
                                c(list(dose=dose),
                                  DLEsamples@data[DLEargNames]))
              
              ## extract the ExpEff function from the model
              EffFun <- slot(Effmodel, "ExpEff")
              
              ## now call the function with dose and with
              ## the arguments taken from the samples
              Effret <- EffFun(dose,Effmodel@data,Effsamples)
              
              ## return the resulting vector
              Gainret <- Effret/(1+(DLEret/(1-DLEret)))
              return(Gainret)
            })


##' @describeIn gain Compute the gain value given a dose level, a pseudo DLE model and a pseudo
##' efficacy model of \code{\linkS4class{Effloglog}} class object without DLE and the efficacy sample
##' @example examples/Model-method-gainNoSamples.R
setMethod("gain",
          signature=
            signature(dose="numeric",
                      DLEmodel="ModelTox",
                      DLEsamples="missing",
                      Effmodel="Effloglog",
                      Effsamples="missing"),
          def=
            function(dose,DLEmodel,Effmodel,...){
              ##extract the prob function from the DLE model
              probFun <- slot(DLEmodel,"prob")
              ##which arguments besides the dose dose it need?
              DLEargNames <- setdiff(names(formals(probFun)),"dose")
              ##now call the function with dose
              DLEvalues<-c()
              for (DLEparName in DLEargNames){
                DLEvalues<-c(DLEvalues, slot(DLEmodel,DLEparName))}
              DLEret <- do.call(probFun,
                                c(list(dose=dose), DLEvalues))
              
              ##extract the ExpEff function from the Eff model
              EffFun <- slot(Effmodel,"ExpEff")
              ##which arguments besides the dose dose it need?
              EffargNames <- setdiff(names(formals(EffFun)),"dose")
              ##now call the function with dose 
              Effvalues<-c()
              for (EffparName in EffargNames){
                Effvalues <- c(Effvalues, slot(Effmodel,EffparName))}
              
              Effret <- do.call(EffFun,
                                c(list(dose=dose),Effvalues))
              Gainret <- Effret/(1+(DLEret/(1-DLEret)))
              return(Gainret)
            })


## ------------------------------------------------------------------------------------
## Update Pseduo models object to obtain new modal estimates for pseudo model parameters
## ------------------------------------------------------------------------------------
## Update the 'LogisticIndepBeta' model
## -----------------------------------------------------------------

##' Update method for the 'LogisticIndepBeta'Model class. This is a method to update the modal
##' estimates of the model parameters \eqn{\phi_1} (phi1) and \eqn{\phi_2} (phi2) when new data 
##' or new observations of responses are available and added in. 
##' 
##' @param object the model of \code{\linkS4class{LogisticIndepBeta}} class object
##' @param data all currently availabvle of \code{\linkS4class{Data}} class object
##' @param \dots unused
##' @return the new \code{\linkS4class{LogisticIndepBeta}} class object
##' 
##' @example examples/Model-method-updateLogisticIndepBeta.R
##' @export
##' @keywords methods
setMethod("update",
          signature=
            signature(object="LogisticIndepBeta"),
          def=
            function(object,
                     data,
                     ...){
              ##Get Pseudo DLE responses (prior) of the model
              
              PseudoDLE<-object@binDLE
              
              ##Get Pseudo DLE weights of the DLE responses of the model
              PseudoDLEweight<-object@DLEweights
              
              
              ##Get the corresponding dose levels for the Pseudo DLE responses from the model
              PseudoDLEdose<- object@DLEdose
              
              ##update the model estimates with data
              model<- LogisticIndepBeta(binDLE=PseudoDLE,DLEweights=PseudoDLEweight,DLEdose=PseudoDLEdose,data=data)
              
              ##return the updated model
              return(model)
            })

## ------------------------------------------------------------------------------------
## Update the 'Effloglog' model
## -----------------------------------------------------------------

##' Update method for the 'Effloglog' Model class. This is a method to update the modal
##' estimates of the model parameters \eqn{\theta_1} (theta1), \eqn{\theta_2} (theta2)  and \eqn{\nu} 
##' (nu, the precision of the efficacy responses) when new data 
##' or new observations of responses are available and added in.
##' 
##' @param object the \code{\linkS4class{Effloglog}} class object
##' @param data all currently available data or responses of \code{\linkS4class{DataDual}}
##' class object
##' @param \dots unused
##' @return the new \code{\linkS4class{Effloglog}} class object
##' 
##' @example examples/Model-method-updateEffloglog.R
##' @export
##' @keywords methods
setMethod("update",
          signature=
            signature(object="Effloglog"),
          def=
            function(object,
                     data,
                     ...){
              ##Get Pseudo Eff responses (prior) of the model
              
              PseudoEff<-object@Eff
              
              ##Get the corresponding dose levels for the Pseudo DLE responses from the model
              PseudoEffdose<- object@Effdose
              
              ## Get the initial values of parameters for nu (if it is not fixed)
              ##OR get the fixed value of nu
              PseudoNu<- object@nu
              
              
              ##update the model estimates with data
              model<- Effloglog(Eff=PseudoEff,Effdose=PseudoEffdose,nu=PseudoNu,data=data)
              
              ##return the updated model
              return(model)
            })
## =================================================================================
## ------------------------------------------------------------------------------------
## Update the 'EffFlexi' model
## -----------------------------------------------------------------

##' Update method for the 'EffFlexi' Model class. This is a method to update 
##' estimates both for the flexible form model and the random walk model (see details in
##' \code{\linkS4class{EffFlexi}} class object) when new data 
##' or new observations of responses are available and added in.
##'  
##' @param object is the model which follow \code{\linkS4class{EffFlexi}} class object
##' @param data all currently available data and responses of \code{\linkS4class{DataDual}}
##' class object
##' @param \dots unused
##' @return the new \code{\linkS4class{EffFlexi}} class object
##' 
##' @example examples/Model-method-updateEffFlexi.R
##' @export
##' @keywords methods
setMethod("update",
          signature=
            signature(object="EffFlexi"),
          def=
            function(object,
                     data,
                     ...){
              ##Get Pseudo Eff responses (prior) of the model
              
              PseudoEff<-object@Eff
              
              ##Get the corresponding dose levels for the Pseudo DLE responses from the model
              PseudoEffdose<- object@Effdose
              
              ## Get the initial values of parameters for Sigma2 (if it is not fixed)
              ##OR get the fixed value of sigma2
              PseudoSigma2<- object@sigma2
              
              
              ## Get the initial values of parameters for Sigma2betaW (if it is not fixed)
              ##OR get the fixed value of sigma2betaW
              PseudoSigma2betaW<- object@sigma2betaW
              
              ##update the model estimates with data
              model<- EffFlexi(Eff=PseudoEff,Effdose=PseudoEffdose,sigma2=PseudoSigma2,sigma2betaW=PseudoSigma2betaW,data=data)
              
              ##return the updated model
              return(model)
            })