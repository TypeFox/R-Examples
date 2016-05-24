#' Interface to MARK for fitting capture-recapture models
#' 
#' Fits user specified models to various types of capture-recapture data by
#' creating input file and running MARK software and retrieving output
#' 
#' This function acts as an interface to the FORTRAN program MARK written by
#' Gary White (\url{http://www.cnr.colostate.edu/~gwhite/mark/mark.htm}). It
#' creates the input file for MARK based on user-specified sub-models
#' (\code{model.parameters}) for each of the parameters in the
#' capture-recapture model being fitted to the data. It runs MARK.EXE (see note
#' below) and then imports the text output file and binary variance-covariance
#' file that were created.  It extracts output values from the text file and
#' creates a list of results that is returned as part of the list (of class
#' mark) which is the return value for this function.
#' 
#' The models that are currently supported are listed in MarkModels.pdf which
#' you can find in the RMark sub-directory of your R Library.  Also, they are
#' listed under Help/Data Types in the MARK interface.
#' 
#' The function mark is a shell that calls 5 other functions in the following
#' order as needed: 1) \code{\link{process.data}}, 2)
#' \code{\link{make.design.data}}, 3) \code{\link{make.mark.model}}, 4)
#' \code{\link{run.mark.model}}, and 5) \code{\link{summary.mark}}. A MARK
#' model can be fitted with this function (\code{mark}) or by calling the
#' individual functions that it uses.  The calling arguments for \code{mark}
#' are a compilation of the calling arguments for each of the functions it
#' calls (with some arguments renamed to avoid conflicts). If data is a
#' processed dataframe (see \code{\link{process.data}}) then it expects to find
#' a value for \code{ddl}.  Likewise, if the data have not been processed, then
#' \code{ddl} should be NULL.  This dual calling structure allows either a
#' single call approach for each model or alternatively for the data to be
#' processed and the design data (\code{ddl}) to be created once and then a
#' whole series of models can be analyzed without repeating those steps.
#' 
#' For descriptions of the arguments \code{data}, \code{begin.time},
#' \code{groups}, \code{age.var}, \code{initial.ages}, \code{age.unit},
#' \code{time.intervals} and \code{mixtures} see \code{\link{process.data}}.
#' 
#' For descriptions of \code{ddl}, \code{design.parameters}=\code{parameters},
#' and \code{right}, see \code{\link{make.design.data}}.
#' 
#' For descriptions of \code{model.name} , \code{model},
#' \code{title},\code{model.parameters}=\code{parameters} ,
#' \code{default.fixed} , \code{initial}, \code{options}, see
#' \code{\link{make.mark.model}}.
#' 
#' And finally, for descriptions of arguments \code{invisible}, \code{filename}
#' and \code{adjust},see \code{\link{run.mark.model}}.
#' 
#' \code{output},\code{silent}, and \code{retry} are the only arguments
#' specific to mark.  \code{output} controls whether a summary of the model
#' input and output are given(if \code{output=TRUE}). \code{silent} controls
#' whether errors are shown when fitting a model. \code{retry} controls the
#' number of times a model will be refitted with new starting values (uses 0)
#' when some parameters are determined to be non-estimable or at a boundary.
#' The latter is the only time it makes sense to retry with new starting values
#' but MARK cannot discern between these two instances. The indices of the beta
#' parameters that are "singular" are stored in \code{results$singular}.
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model.name Optional name for the model
#' @param model Type of c-r model (eg CJS, Burnham, Barker)
#' @param title Optional title for the MARK analysis output
#' @param model.parameters List of model parameter specifications
#' @param initial Optional vector of named or unnamed initial values for beta
#' parameters or previously run model object
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param right if TRUE, any intervals created in design.parameters are closed
#' on the right and open on left and vice-versa if FALSE
#' @param groups Vector of names factor variables for creating groups
#' @param age.var Optional index in groups vector of a variable that represents
#' age
#' @param initial.ages Optional vector of initial ages for each age level
#' @param age.unit Increment of age for each increment of time
#' @param time.intervals Intervals of time between the capture occasions
#' @param nocc number of occasions for Nest model; either time.intervals or
#' nocc must be specified for this model
#' @param output If TRUE produces summary of model input and model output
#' @param invisible If TRUE, window for running MARK is hidden
#' @param adjust If TRUE, adjusts npar to number of cols in design matrix,
#' modifies AIC and records both
#' @param mixtures number of mixtures for heterogeneity model or number of secondary samples for MultScaleOcc model
#' @param se if TRUE, se and confidence intervals are shown in summary sent to
#' screen
#' @param filename base filename for files created by MARK.EXE. Files are named
#' filename.*.
#' @param prefix base filename prefix for files created by MARK.EXE; for
#' example if prefix="SpeciesZ" files are named "SpeciesZnnn.*"
#' @param default.fixed if TRUE, real parameters for which the design data have
#' been deleted are fixed to default values
#' @param silent if TRUE, errors that are encountered are suppressed
#' @param retry number of reanalyses to perform with new starting values when
#' one or more parameters are singular
#' @param options character string of options for Proc Estimate statement in
#' MARK .inp file
#' @param brief if TRUE and output=TRUE then a brief summary line is given
#' instead of a full summary for the model
#' @param realvcv if TRUE the vcv matrix of the real parameters is extracted
#' and stored in the model results
#' @param delete if TRUE the output files are deleted after the results are
#' extracted
#' @param external if TRUE the mark object is saved externally rather than in
#' the workspace; the filename is kept in its place
#' @param profile.int if TRUE will compute profile intervals for each real
#' parameter; or you can specify a vector of real parameter indices
#' @param chat value of chat used for profile intervals
#' @param reverse if set to TRUE, will reverse timing of transition (Psi) and
#' survival (S) in Multistratum models
#' @param run if FALSE does not run model after creation 
#' @param input.links specifies set of link functions for parameters with non-simplified structure
#' @param parm.specific if TRUE, forces a link to be specified for each parameter
#' @param mlogit0 if TRUE, any real parameter that is fixed to 0 and has an mlogit link will 
#' have its link changed to logit so it can be simplified
#' @param threads number of cpus to use with mark.exe if positive or number of cpus to remain idle if negative
#' @param hessian if TRUE specifies to MARK to use hessian rather than second partial matrix
#' @param accumulate if TRUE accumulate like data values into frequencies
#' @param allgroups Logical variable; if TRUE, all groups are created from
#' factors defined in \code{groups} even if there are no observations in the
#' group
#' @param strata.labels vector of single character values used in capture
#' history(ch) for ORDMS, CRDMS, RDMSOccRepro models; it can contain one more value beyond what is
#' in ch for an unobservable state except for RDMSOccRepro which is used to specify strata ordering (eg 0 not-occupied, 1 occupied no repro, 2 occupied with repro.
#' @param counts named list of numeric vectors (one group) or matrices (>1
#' group) containing counts for mark-resight models
#' @param icvalues numeric vector of individual covariate values for computation of real values
#' @param wrap if TRUE, data lines are wrapped to be length 80; if length of a row is not a 
#'   problem set to FALSE and it will run faster
#' @return model: a MARK object containing output and extracted results. It is
#' a list with the following elements \item{data}{name of the processed data
#' frame} \item{model}{type of analysis model (see list above)}
#' \item{title}{title used for analysis} \item{model.name}{descriptive name of
#' model} \item{links}{vector of link function(s) used for parameters, one for
#' each row in design matrix or only one if all parameters use the same
#' function} \item{mixtures}{number of mixtures in Pledger-style closed
#' capture-recapture models} \item{call}{call to make.mark.model used to
#' construct the model} \item{parameters}{a list of parameter descriptions
#' including the formula, pim.type, link etc.} \item{model.parameters}{the list
#' of parameter descriptions used in the call to mark; this is used only by
#' \code{rerun.mark}} \item{time.intervals}{Intervals of time between the
#' capture occasions} \item{number.of.groups}{number of groups defined in the
#' data} \item{group.labels}{vector of labels for the groups}
#' \item{nocc}{number of capture occasions} \item{begin.time}{single time of
#' vector of times (if different for groups) for the first capture occasion}
#' \item{covariates}{vector of covariate names (as strings) used in the model}
#' \item{fixed}{dataframe of parameters set at fixed values; \code{index} is
#' the parameter index in the full parameter structure and \code{value} is the
#' fixed value for the real parameter} \item{design.matrix}{design matrix used
#' in the input to MARK.EXE} \item{pims}{list of pims used for each parameter
#' including any group or strata designations; each parameter in list is
#' denoted by name and within each parameter one or more sub-lists represent
#' groups and strata if any} \item{design.data}{design data used to construct
#' the design matrix} \item{strata.labels}{labels for strata if any}
#' \item{mlogit.list}{structure used to simplify parameters that use mlogit
#' links} \item{simplify}{list containing \code{pim.translation} which
#' translate between all different and simplified pims, \code{real.labels}
#' which are labels for real parameters for full (non-simplified) pim structure
#' and \code{links} the link function names for the full parameter structure}
#' \item{output}{base portion of filenames for input,output, vc and residual
#' files output from MARK.EXE} \item{results}{List of values extracted from
#' MARK ouput} \tabular{lll}{ \tab \code{lnl} \tab -2xLog Likelihood value \cr
#' \tab \code{npar} \tab Number of parameters (always the number of columns in
#' design matrix) \cr \tab \code{npar.unadjusted} \tab number of estimated
#' parameters from MARK if different than npar \cr \tab \code{n} \tab effective
#' sample size \cr \tab \code{AICc} \tab Small sample corrected AIC using npar
#' \cr \tab \code{AICc.unadjusted} \tab Small sample corrected AIC using
#' npar.unadjusted \cr \tab \code{beta} \tab data frame of beta parameters with
#' estimate, se, lcl, ucl \cr \tab \code{real} \tab data frame of real
#' parameters with estimate, se, lcl, ucl and fixed \cr \tab \code{beta.vcv}
#' \tab variance-covariance matrix for beta \cr \tab \code{derived} \tab
#' dataframe of derived parameters if any \cr \tab \code{derived.vcv} \tab
#' variance-covariance matrix for derived parameters if any \cr \tab
#' \code{covariate.values} \tab dataframe with fields \code{Variable} and
#' \code{Value} \cr \tab \tab which are the covariate names and value used for
#' real parameter \cr \tab \tab estimates in the MARK output\cr \tab
#' \code{singular} \tab indices of beta parameters that are non-estimable or at
#' a boundary \cr \tab \code{real.vcv} \tab variance-covariance matrix for real
#' parameters (simplified) if realvcv=TRUE \cr } \item{chat}{over-dispersion
#' constant; if not present assumed to be 1}
#' @note It is assumed that MARK.EXE is located in directory "C:/Program
#' Files/Mark".  If it is in a different location set the variable MarkPath to
#' the directory location. For example, seting MarkPath="C:/Mark/" at the R
#' prompt will assign run "c:/mark/mark.exe" to do the analysis.  If you have
#' chosen a non-default path for Mark.exe, MarkPath needs to be defined for
#' each R session.  It is easiest to do this assignment automatically by
#' putting the MarkPath assignment into your .First function which is run each
#' time an R session is initiated.  In addition to MarkPath, the variable
#' MarkViewer can be assigned to a program other than notepad.exe (see
#' \code{\link{print.mark}}).
#' @author Jeff Laake
#' @export
#' @importFrom stats as.formula coef formula median model.matrix optim 
#'            optimize plogis pnorm pt qchisq qnorm terms uniroot
#' @importFrom utils read.delim read.fwf read.table
#'              type.convert write.table
#' @seealso \code{\link{make.mark.model}}, \code{\link{run.mark.model}},
#' \code{\link{make.design.data}}, \code{\link{process.data}},
#' \code{\link{summary.mark}}
#' @keywords models
#' @examples
#' 
#' data(dipper)
#' dipper.Phidot.pdot=mark(dipper,threads=1)
#' 
mark <-
function(data,ddl=NULL,begin.time=1,model.name=NULL,model="CJS",title="",model.parameters=list(),initial=NULL,
design.parameters=list(), right=TRUE, groups = NULL, age.var = NULL, initial.ages = 0, age.unit = 1, time.intervals = NULL,nocc=NULL,output=TRUE,
invisible=TRUE,adjust=TRUE,mixtures=1,se=FALSE,filename=NULL,prefix="mark",default.fixed=TRUE,silent=FALSE,retry=0,options=NULL,brief=FALSE,
realvcv=FALSE,delete=FALSE,external=FALSE,profile.int=FALSE,chat=NULL,reverse=FALSE,run=TRUE,input.links=NULL,parm.specific=FALSE,mlogit0=FALSE,threads=-1,hessian=FALSE,accumulate=TRUE,
allgroups=FALSE,strata.labels=NULL,counts=NULL,icvalues=NULL,wrap=TRUE)
{
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
# 
simplify=TRUE
if(is.null(data$data))
{
   if(!is.null(ddl))
   {
      message("Warning: specification of ddl ignored, as data have not been processed\n")
      ddl=NULL
   }
   data.proc=process.data(data,begin.time=begin.time, model=model,mixtures=mixtures, 
                          groups = groups, age.var = age.var, initial.ages = initial.ages, 
                          age.unit = age.unit, time.intervals = time.intervals,nocc=nocc,reverse=reverse,
				          allgroups=allgroups, strata.labels=strata.labels,counts=counts)
}   
else
   data.proc=data
#
# If the design data have not been constructed, do so now
#
if(is.null(ddl)) ddl=make.design.data(data.proc,design.parameters,right=right)
#
#  check to make sure all entered as lists
#
tryCatch(length(model.parameters), error = function(e) message("Make sure you have a tilde at the beginning of each formula\n"))
if(length(model.parameters)!=0)
	for(i in 1:length(model.parameters))
	{
		if(!is.list(model.parameters[[i]]))
			stop("\nEach parameter distribution must be specified as a list\n")
		if(is.language(model.parameters[[i]][[1]])&(is.null(names(model.parameters[[i]])) || names(model.parameters[[i]])[1]==""))
			message("Make sure you have an = between formula and tilde for formula\n")
	}
#
# Run model as many as times requested if needed
#
i=0
converge=FALSE
while(i<=retry & !converge)
{
#
# Make the model with specified or default parameters
#
   if(is.list(model.parameters))
   {
      model<-try(make.mark.model(data.proc,title=title,parameters=model.parameters,
             ddl=ddl,initial=initial,call=match.call(),default.fixed=default.fixed,
             model.name=model.name,options=options,profile.int=profile.int,chat=chat,
			 input.links=input.links,parm.specific=parm.specific,mlogit0=mlogit0,hessian=hessian,
			 accumulate=accumulate,icvalues=icvalues,wrap=wrap))
      if(class(model)[1]=="try-error")
	  {
		  stop("Misspecification of model or internal error in code")
	  }
	  else
         model$model.parameters=model.parameters
	  if(!run)return(model)
   }
   else
      stop("Model parameters must be specified as a list")
#
# Summarize model input if output=TRUE
#
   if(output & i==1)
   {
     cat("\n")
     print(summary(model))
   }
#
# Run model
#
   if(silent)
	   runmodel<-suppressMessages(run.mark.model(model,invisible=invisible,adjust=adjust,filename=filename,prefix=prefix,realvcv=realvcv,delete=delete,threads=threads,ignore.stderr=silent))
   else
       runmodel<-run.mark.model(model,invisible=invisible,adjust=adjust,filename=filename,prefix=prefix,realvcv=realvcv,delete=delete,threads=threads,ignore.stderr=silent)
   if(is.null(runmodel))
   {
     if(!silent)message("\n\n********Following model failed to run :",model$model.name,"**********\n\n")
     return(invisible())
   }
   else
   {
#
#  Check if any parameters are singular and if retry=TRUE, the refit model with
#  new initial values
#
      if(retry>0 && !is.null(runmodel$results$singular))
      {
		 if(!silent)message("\nRe-running analysis with new starting values\n")
         i=i+1
         converge=FALSE
         initial=runmodel$results$beta$estimate
		 names(initial)=rownames(runmodel$results$beta)
         initial[runmodel$results$singular]=0
         next
      }
      else
         converge=TRUE
    }
}
#
# Summarize model results if output=TRUE
#
   if(output)
   {
      if(!brief)
      {
         cat("\n")
         print(summary(runmodel,se=se))
      }
      else
      {
         sum.res=summary(runmodel)
         cat(paste("\n Model:",sum.res$model.name," npar=",sum.res$npar," lnl = ",sum.res$lnl,"AICc =",sum.res$AICc))
      }
   }
#
# Return fitted MARK model object or if external, return character string with same class and save file
#
   if(external)
   {
     marksave=paste(runmodel$output,".rda",sep="")
     model=runmodel
     save(model,file=marksave)
     class(marksave)=class(runmodel)
     return(marksave)
   } else
     return(runmodel)
}
