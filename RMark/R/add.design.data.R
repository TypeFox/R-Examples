#' Add design data
#' 
#' Creates new design data fields in (\code{ddl}) that bin the fields
#' \code{cohort}, \code{age} or \code{time}. Other fields (e.g., effort value
#' for time) can be added to \code{ddl} with R commands.
#' 
#' Design data can be added to the parameter specific design dataframes with R
#' commands.  Often the additional fields will be functions of \code{cohort},
#' \code{age} or \code{time}.  \code{add.design.data} provides an easy way to
#' add fields that bin (put into intervals) the original values of
#' \code{cohort}, \code{age} or \code{time}.  For example, \code{age} may have
#' levels from 0 to 10 which means the formula \code{~age} will have 11
#' parameters, one for each level of the factor.  It might be more desirable
#' and more parimonious to have a simpler 2 age class model of young and
#' adults.  This can be done easily by adding a new design data field that bins
#' \code{age} into 2 intervals (age 0 and 1+) as in the following example:
#' 
#' \preformatted{ ddl=make.design.data(proc.example.data)
#' ddl=add.design.data(proc.example.data,ddl,parameter="Phi",type="age",
#' bins=c(0,.5,10),name="2ages") }
#' 
#' By default, the bins are open on the left and closed on the right (i.e.,
#' binning x by (x1,x2] is equivalent to x1<x<=x2) except for the first
#' interval which is closed on the left.  Thus, for the above example, the age
#' bins are [0,.5] and (.5,10].  Since the ages in the example are 0,1,2...
#' using any value >0 and <1 in place of 0.5 would bin the ages into 2 classes
#' of 0 and 1+.  This behavior can be modified by changing the argument
#' right=FALSE to create an interval that is closed on the left and open on the
#' right.  In some cases this can make reading the values of the levels
#' somewhat easier. It is important to recognize that the new variable is only
#' added to the design data for the defined \code{parameter} and can only be
#' used in model formula for that parameter.  Multiple calls to
#' \code{add.design.data} can be used to add the same or different fields for
#' the various parameters in the model. For example, the same 2 age class
#' variable can be added to the design data for p with the command:
#' 
#' \preformatted{
#' ddl=add.design.data(proc.example.data,ddl,parameter="p",type="age",
#' bins=c(0,.5,10),name="2ages") }
#' 
#' The \code{name} must be unique within the parameter design data, so they
#' should not use pre-defined values of \code{group, age, Age, time, Time,
#' cohort, Cohort}.  If you choose a \code{name} that already exists in the
#' design data for the \code{parameter}, it will not be added but it can
#' replace the variable if \code{replace=TRUE}. For example, the \code{2ages}
#' variable can be re-defined to use 0-1 and 2+ with the command:
#' 
#' \preformatted{
#' ddl=add.design.data(proc.example.data,ddl,parameter="Phi",type="age",
#' bins=c(0,1,10),name="2ages",replace=TRUE) }
#' 
#' Keep in mind that design data are stored with the \code{mark} model object
#' so if a variable is redefined, as above, this could become confusing if some
#' models have already been constructed using a different definition for the
#' variable. The model formula and names would appear to be identical but they
#' would have a different model structure.  The difference would be apparent if
#' you examined the design data and design matrix of the model object but would
#' the difference would be transparent based on the model names and formula.
#' Thus, it would be best to avoid constructing models from design data fields
#' with different structures but the same name.
#' 
#' @param data processed data list resulting from \code{\link{process.data}}
#' @param ddl current design dataframe initially created with
#' \code{\link{make.design.data}}
#' @param parameter name of model parameter (e.g., "Phi" for CJS models)
#' @param type either "age", "time" or "cohort"
#' @param bins bins for grouping
#' @param name name assigned to variable in design data
#' @param replace if TRUE, replace any variable with same name as \code{name}
#' @param right If TRUE, bin intervals are closed on the right
#' @export
#' @return Design data list with new field added for the specified parameter.
#' See \code{\link{make.design.data}} for a description of the list structure.
#' @note For the specific case of "closed" capture models, the parameters
#' \code{p} (capture probability) and \code{c} (recapture probability) can be
#' treated in a special fashion.  Because they really the same type of
#' parameter, it is useful to be able to share a common model structure (i.e.,
#' same columns in the design matrix). This is indicated with the
#' \code{share=TRUE} element in the model description for \code{p}.  If the
#' parameters are shared then the additional covariate \code{c} is added to the
#' design data, which is \code{c=0} for parameter \code{p} and \code{c=1} for
#' parameter \code{c}.  This enables an additive model to be developed where
#' recapture probabilities mimic the pattern in capture probabilities except
#' for an additive constant.  The covariate \code{c} can only be used in the
#' model for \code{p} if \code{share=TRUE}. If the latter is not set using
#' \code{c} in a formula will result in an error.  Likewise, if
#' \code{share=TRUE}, then the design data for \code{p} and \code{c} must be
#' the same because the design data are merged in constructing the design
#' matrix.  Thus if you add design data for parameter \code{p}, you should add
#' a similar field for parameter \code{c} if you intend to fit shared models
#' for the two parameters.  If the design data do not match and you try to fit
#' a shared model, an error will result.
#' @author Jeff Laake
#' @seealso \code{\link{make.design.data}}, \code{\link{process.data}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(example.data)
#' example.data.proc=process.data(example.data)
#' ddl=make.design.data(example.data.proc)
#' ddl=add.design.data(example.data.proc,ddl,parameter="Phi",type="age",
#'   bins=c(0,.5,10),name="2ages")
#' ddl=add.design.data(example.data.proc,ddl,parameter="p",type="age",
#' bins=c(0,.5,10),name="2ages")
#' ddl=add.design.data(example.data.proc,ddl,parameter="Phi",type="age",
#' bins=c(0,1,10),name="2ages",replace=TRUE)
#' }
add.design.data <-
function(data,ddl,parameter,type="age",bins=NULL,name=NULL,replace=FALSE,right=TRUE)
# -------------------------------------------------------------------------------------------------------------
#
# add.design.data  - enables design data fields to be added to the design.data(ddl)
#                    the added fields must be a function of cohort, age or time.
#                    Any other fields (eg effort value for time) can be added with R
#                    commands, but this may be changed later. 
#
# Arguments:
#
# data             - data list resulting from process.data 
# ddl              - current design dataframe 
# parameter        - parameter name
# type             - either "age", "time" or "cohort"
# bins             - bins for grouping 
# name             - variable name
# replace	   - if TRUE, replace field with same name
#
# Value:  
#
#  ddl - modified design data 
#
#
# Functions used: compute.design.data, setup.parameters, valid.parameters
#
# -------------------------------------------------------------------------------------------------------------
{
#
# Check validity of parameter list, if any given
#
  if(!valid.parameters(data$model,parameter))stop()
#
# Make sure that a name has been given for added variable
#
  if(is.null(name)) stop("A name is required for an added design variable")
#
#  Setup variables in parameter list
#
  parameters=setup.parameters(data$model)
  model.list=setup.model(data$model,data$nocc,data$mixtures)
#
#  Assign user defined pim.types
#
  for(pname in names(ddl$pimtypes))
     parameters[[pname]]$pim.type=ddl$pimtypes[[pname]]$pim.type
#
#  Compute design data for the parameter
#
  if(!model.list$robust) parameters[[parameter]]$secondary=FALSE
#
# Compute design data for this parameter
#
  if(data$mixtures==1)
  {
     parameters$mix=FALSE
     parameters[[parameter]]$rows=1
  }
  if(!is.null(parameters[[parameter]]$bystratum) && parameters[[parameter]]$bystratum)
  {
     strata.labels=data$strata.labels
     nstrata=data$nstrata
     if(!is.null(parameters[[parameter]]$tostrata) && parameters[[parameter]]$tostrata)
        tostrata=TRUE
     else
        tostrata=FALSE
  }
  else
  {
     strata.labels=NULL
     nstrata=1
     tostrata=FALSE
  }
#
# Compute design data
#
  design.data=compute.design.data(data,parameters[[parameter]]$begin,parameters[[parameter]]$num,
                   parameters[[parameter]]$type,parameters[[parameter]]$mix,parameters[[parameter]]$rows,
                   parameters[[parameter]]$pim.type,parameters[[parameter]]$secondary, nstrata,
                   tostrata,strata.labels)
#
#  Limit design data to rows in ddl to handle elimination of design data
#
   design.data=design.data[row.names(design.data)%in%row.names(ddl[[parameter]]),]
#
# Add variable depending on type
#

  if(type=="cohort")    
  {
     if(is.null(bins))
        new.data=as.factor(design.data$cohort)
     else
        new.data=cut(design.data$cohort,bins,include.lowest=TRUE,right=right)
  } else
  if(type=="age")
  {
     if(is.null(bins))
        new.data=as.factor(design.data$age)
     else
        new.data=cut(design.data$age,bins,include.lowest=TRUE,right=right)
  } else
  if(type=="time")
  {
     if(is.null(bins)) 
        new.data=as.factor(design.data$time)
     else
        new.data=cut(design.data$time,bins,include.lowest=TRUE,right=right)
  } else stop("invalid type")
  vnames=names(ddl[[parameter]])
  if(name %in% vnames)
  {
    if(replace)
       ddl[[parameter]][,name]=new.data 
    else
       stop(paste("Variable ",name," already in design data. Use replace=TRUE if you want to replace current values"))
  }
  else
  {
     ddl[[parameter]]<-cbind(ddl[[parameter]],new.data)
     names(ddl[[parameter]])<- c(vnames,name)
  }
return(ddl)
}
