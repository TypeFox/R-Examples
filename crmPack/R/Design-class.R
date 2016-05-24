#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com]
##         Wai Yin Yeung [ w*.* yeung1 *a*t* lancaster *.* ac *.* uk]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Design-class.R] by DSB Son 18/01/2015 21:35>
##
## Description:
## This class encapsulates a whole CRM design.
##
## History:
## 12/02/2014   file creation
## 10/07/2015  adding designs for Pseudo models
#####################################################################################

##' @include Model-class.R
##' @include Rules-class.R
##' @include Data-class.R
##' @include helpers.R
{}


## --------------------------------------------------
## Classes for rule-based designs
## --------------------------------------------------


##' Class for rule-based designs
##'
##' The difference to \code{\linkS4class{Design}} class is that
##' model, stopping and increments slots are missing.
##'
##' @slot nextBest how to find the next best dose, an object of class
##' \code{\linkS4class{NextBest}}
##' @slot cohortSize rules for the cohort sizes,
##' an object of class \code{\linkS4class{CohortSize}}
##' @slot data what is the dose grid, any previous data, etc., contained
##' in an object of class \code{\linkS4class{Data}}
##' @slot startingDose what is the starting dose? Must lie on the grid in
##' \code{data}
##'
##' @example examples/design-class-RuleDesign.R
##' @export
##' @keywords classes
.RuleDesign <-
    setClass(Class="RuleDesign",
             representation(nextBest="NextBest",
                            cohortSize="CohortSize",
                            data="Data",
                            startingDose="numeric"),
             prototype(nextBest=.NextBestThreePlusThree(),
                       cohortSize=CohortSizeConst(3),
                       data=Data(doseGrid=1:3),
                       startingDose=1),
             validity=
                 function(object){
                     o <- Validate()

                     o$check(is.scalar(object@startingDose),
                             "startingDose must be scalar")
                     o$check(object@startingDose %in% object@data@doseGrid,
                             "startingDose must be included in data@doseGrid")

                     o$result()
                 })
validObject(.RuleDesign())

##' Initialization function for "RuleDesign"
##'
##' @param nextBest see \code{\linkS4class{RuleDesign}}
##' @param cohortSize see \code{\linkS4class{RuleDesign}}
##' @param data see \code{\linkS4class{RuleDesign}}
##' @param startingDose see \code{\linkS4class{RuleDesign}}
##' @return the \code{\linkS4class{RuleDesign}} object
##'
##' @export
##' @keywords methods
RuleDesign <- function(nextBest,
                       cohortSize,
                       data,
                       startingDose)
{
    .RuleDesign(nextBest=nextBest,
                cohortSize=cohortSize,
                data=data,
                startingDose=as.numeric(startingDose))
}


## --------------------------------------------------
## Classes for model-based designs
## --------------------------------------------------

##' Class for the CRM design
##'
##' In addition to the slots in the more simple \code{\linkS4class{RuleDesign}},
##' objects of this class contain:
##'
##' @slot model the model to be used, an object of class
##' \code{\linkS4class{Model}}
##' @slot stopping stopping rule(s) for the trial, an object of class
##' \code{\linkS4class{Stopping}}
##' @slot increments how to control increments between dose levels,
##' an object of class \code{\linkS4class{Increments}}
##' @slot PLcohortSize rules for the cohort sizes for placebo, if any planned
##' an object of class \code{\linkS4class{CohortSize}}
##'
##' @example examples/design-class-Design.R
##' @export
##' @keywords classes
.Design <-
    setClass(Class="Design",
             representation(model="Model",
                            stopping="Stopping",
                            increments="Increments",
                            PLcohortSize="CohortSize"),
             prototype(model=.LogisticNormal(),
                       nextBest=.NextBestNCRM(),
                       stopping=.StoppingMinPatients(),
                       increments=.IncrementsRelative(),
                       PLcohortSize=CohortSizeConst(1)),
             contains=list("RuleDesign"))
validObject(.Design())


##' Initialization function for "Design"
##'
##' @param model see \code{\linkS4class{Design}}
##' @param stopping see \code{\linkS4class{Design}}
##' @param increments see \code{\linkS4class{Design}}
##' @param PLcohortSize see \code{\linkS4class{Design}}
##' @param \dots additional arguments for \code{\link{RuleDesign}}
##' @return the \code{\linkS4class{Design}} object
##'
##' @export
##' @keywords methods
Design <- function(model,
                   stopping,
                   increments,
                   PLcohortSize=CohortSizeConst(1),
                   ...)
{
    start <- RuleDesign(...)
    .Design(start,
            model=model,
            stopping=stopping,
            increments=increments,
            PLcohortSize=PLcohortSize)
}



##' Class for the dual-endpoint CRM design
##'
##' This class has special requirements for the \code{model} and \code{data}
##' slots in comparison to the parent class \code{\linkS4class{Design}}:
##'
##' @slot model the model to be used, an object of class
##' \code{\linkS4class{DualEndpoint}}
##' @slot data what is the dose grid, any previous data, etc., contained
##' in an object of class \code{\linkS4class{DataDual}}
##'
##' Note that the \code{NextBest} slot can be of any class, this allows for easy
##' comparison with recommendation methods that don't use the
##' biomarker information.
##'
##' @example examples/design-class-DualDesign.R
##' @export
##' @keywords classes
.DualDesign <-
    setClass(Class="DualDesign",
             representation(model="DualEndpoint",
                            data="DataDual"),
             prototype(model=.DualEndpoint(),
                       nextBest=.NextBestDualEndpoint(),
                       data=DataDual(doseGrid=1:2),
                       startingDose=1),
             contains=list("Design"))
validObject(.DualDesign())


##' Initialization function for "DualDesign"
##'
##' @param model see \code{\linkS4class{DualDesign}}
##' @param data see \code{\linkS4class{DualDesign}}
##' @param \dots additional arguments for \code{\link{Design}}
##' @return the \code{\linkS4class{DualDesign}} object
##'
##' @export
##' @keywords methods
DualDesign <- function(model,
                       data,
                       ...)
{
    start <- Design(data=data,
                    model=model,
                    ...)
    .DualDesign(start,
                model=model,
                data=data)
}




##' Creates a new 3+3 design object from a dose grid
##'
##' @param doseGrid the dose grid to be used
##' @return the object of class \code{\linkS4class{RuleDesign}} with the
##' 3+3 design
##'
##' @example examples/design-class-ThreePlusThreeDesign.R
##' @export
##' @keywords programming
##' @author Daniel Sabanes Bove \email{sabanesd@@roche.com}
ThreePlusThreeDesign <- function(doseGrid)
{
    emptydata <- Data(doseGrid=doseGrid)

    design <- RuleDesign(nextBest=NextBestThreePlusThree(),
                         data=emptydata,
                         cohortSize=CohortSizeConst(size=3L),
                         ## using a constant cohort size of 3,
                         ## we obtain exactly the 3+3 design
                         startingDose=head(emptydata@doseGrid, 1))

    return(design)
}

## ===================================================================================
## -------------------------------------------------------------------------------
## Design class using DLE responses only based on the pseudo DLE model with samples
## ---------------------------------------------------------------------------
##' This is a class of design based only on DLE responses using the 'LogisticIndepBeta' class model
##' and DLE samples are also used. 
##' In addition to the slots in the more simple \code{\linkS4class{RuleDesign}},
##' objects of this class contain:
##' 
##' @slot model the pseudo DLE model to be used, an object class of 
##' \code{\linkS4class{ModelTox}}
##' @slot stopping stopping rule(s) for the trial, an object class of \code{\linkS4class{Stopping}}
##' @slot increments how to control increments between dose levels, an object class of 
##' \code{\linkS4class{Increments}}
##' 
##' @example examples/design-class-TDsamplesDesign.R
##' @export
##' @keywords class
.TDsamplesDesign <-
  setClass(Class="TDsamplesDesign",
           representation(model="ModelTox",
                          stopping="Stopping",
                          increments="Increments"),
           prototype(model=.LogisticIndepBeta(),
                     nextBest=.NextBestTDsamples(),
                     stopping=.StoppingMinPatients(),
                     increments=.IncrementsRelative()),
           contains=list("RuleDesign"))

validObject(.TDsamplesDesign())
##' Initialization function for 'TDsamplesDesign' class
##' 
##' @param model see \code{\linkS4class{TDsamplesDesign}}
##' @param stopping see \code{\linkS4class{TDsamplesDesign}}
##' @param increments see \code{\linkS4class{TDsamplesDesign}}
##' @param \dots additional arguments for \code{\linkS4class{RuleDesign}}
##' @return the \code{\linkS4class{TDsamplesDesign}} class object
##' 
##' @export
##' @keywords methods
TDsamplesDesign<-function(model,stopping,increments,...){
  start<-RuleDesign(...)
  .TDsamplesDesign(start,model=model,stopping=stopping,increments=increments)
}

## =============================================================================
## -------------------------------------------------------------------------------
##' Design class using DLE responses only based on the pseudo DLE model without sample
##'
##' This is a class of design based only on DLE responses using the 'LogisticIndepBeta' class model
##' are used without samples.
##' In addition to the slots in the more simple \code{\linkS4class{RuleDesign}},
##' objects of this class contain:
##' 
##' @slot model the pseudo DLE model to be used, an object class of 
##' \code{\linkS4class{ModelTox}}
##' @slot stopping stopping rule(s) for the trial, an object class of \code{\linkS4class{Stopping}}
##' @slot increments how to control increments between dose levels, an object class of 
##' \code{\linkS4class{Increments}}
##' 
##' @example examples/design-class-TDDesign.R
##' @export
##' @keywords class 
.TDDesign <-
  setClass(Class="TDDesign",
           representation(model="ModelTox",
                          stopping="Stopping",
                          increments="Increments"),
           prototype(model=.LogisticIndepBeta(),
                     nextBest=.NextBestTD(),
                     stopping=.StoppingMinPatients(),
                     increments=.IncrementsRelative()),
           contains=list("RuleDesign"))

validObject(.TDDesign())
##' Initialization function for 'TDDesign' class
##' 
##' @param model please refer to \code{\linkS4class{TDDesign}} class object
##' @param stopping please refer to \code{\linkS4class{TDDesign}} class object
##' @param increments please refer to \code{\linkS4class{TDDesign}} class object
##' @param \dots additional arguments for \code{\linkS4class{RuleDesign}}
##' @return the \code{\linkS4class{TDDesign}} class object
##' 
##' @export
##' @keywords methods
TDDesign<-function(model,stopping,increments,...){
  start<-RuleDesign(...)
  .TDDesign(start,model=model,stopping=stopping,increments=increments)}



## ---------------------------------------------------------------------------------------------------
## class for design based on DLE and efficacy response with samples using pseudo DLE and efficacy models
##----------------------------------------------------------------------------------------------------
##' This is a class of design based on DLE responses using the \code{\linkS4class{LogisticIndepBeta}} model 
##' model and efficacy responses using \code{\linkS4class{ModelEff}}  model class
##' with DLE and efficacy samples.It contain all slots in 
##' \code{\linkS4class{RuleDesign}} and \code{\linkS4class{TDsamplesDesign}} class object
##' 
##' @slot data the data set of \code{\linkS4class{DataDual}} class object
##' @slot Effmodel the pseudo efficacy model to be used, an object class of 
##' \code{\linkS4class{ModelEff}}
##' 
##' @example examples/design-class-DualResponsesSamplesDesign.R
##' @export
##' @keywords class 
##' 
.DualResponsesSamplesDesign <-
setClass(Class="DualResponsesSamplesDesign",
         representation(Effmodel="ModelEff",
                        data="DataDual"),
         prototype(nextBest=.NextBestMaxGainSamples(),
                   data=DataDual(doseGrid=1:2),
                   startingDose=1,
                   model=.LogisticIndepBeta()),
         contains=list("TDsamplesDesign")
           )
validObject(.DualResponsesSamplesDesign())

##' Initialization function for 'DualResponsesSamplesDesign"
##' @param data please refer to \code{\linkS4class{DualResponsesSamplesDesign}} class object
##' @param Effmodel please refer to \code{\linkS4class{DualResponsesSamplesDesign}} class object
##' @param \dots additional arguments for \code{\link{TDsamplesDesign}}
##' 
##' @return the \code{\linkS4class{DualResponsesSamplesDesign}} class object
##' 
##' @export
##' @keywords methods
DualResponsesSamplesDesign <- function(Effmodel,
                                data,
                                ...)
  
{
  
  start <- TDsamplesDesign(data=data,...)
  .DualResponsesSamplesDesign(start,
                              Effmodel=Effmodel,
                              data=data)
}

## ---------------------------------------------------------------------------------------------------
## class for design based on DLE and efficacy response without  samples using pseudo DLE and efficacy models
##----------------------------------------------------------------------------------------------------
##' This is a class of design based on DLE responses using the \code{\linkS4class{LogisticIndepBeta}} model 
##' model and efficacy responses using \code{\linkS4class{ModelEff}}  model class
##' without DLE and efficacy samples. It contain all slots in 
##' \code{\linkS4class{RuleDesign}} and \code{\linkS4class{TDDesign}} class object
##' 
##' @slot data the data set of \code{\linkS4class{DataDual}} class object
##' @slot Effmodel the pseudo efficacy model to be used, an object class of 
##' \code{\linkS4class{ModelEff}}
##' 
##' @example examples/design-class-DualResponsesDesign.R
##' @export
##' @keywords class 
.DualResponsesDesign <-
  setClass(Class="DualResponsesDesign",
           representation(Effmodel="ModelEff",
                          data="DataDual"),
           prototype(nextBest=.NextBestMaxGain(),
                     data=DataDual(doseGrid=1:2),
                     startingDose=1,
                     model=.LogisticIndepBeta()),
           contains=list("TDDesign")
  )
validObject(.DualResponsesDesign())


##' Initialization function for 'DualResponsesDesign"
##' @param data please refer to \code{\linkS4class{DualResponsesDesign}} class object
##' @param Effmodel please refer to \code{\linkS4class{DualResponsesDesign}} class object
##' @param \dots additional arguments for \code{\link{TDDesign}}
##' @return the \code{\linkS4class{DualResponsesDesign}} class object
##' 
##' @export
##' @keywords methods
DualResponsesDesign <- function(Effmodel,
                                data,
                                ...)
  
{
  
  start <- TDDesign(data=data,...)
  .DualResponsesDesign(start,
                       Effmodel=Effmodel,
                       data=data)
}

  ## ===============================================================================