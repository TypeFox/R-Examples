#' Runs a previous MARK model with new starting values
#' 
#' Runs a previous MARK model with new starting values but without specifying
#' the model parameter formulas.  This function is most useful with
#' \code{mark.wrapper} in which a list of models is analyzed and the set of
#' formulas are not specified for each model.
#' 
#' This is a simple function that restarts an analysis with MARK typically
#' using another model for initial values of the beta parameters.  The
#' processed dataframe (\code{data}) and design data list (\code{ddl}) must be
#' specified but the \code{model.parameters} are extracted from \code{model}.
#' \code{initial} values are not optional otherwise this would be no different
#' than the original call to \code{mark}. More complete definitions of the
#' arguments can be found in \code{\link{mark}} or \code{\link{run.mark.model}}
#' or \code{\link{make.mark.model}}.
#' 
#' @param model previously run MARK model
#' @param data processed dataframe used with model
#' @param ddl design data list used with model
#' @param initial vector of initial values for beta parameters or previously
#' run model object of similar structure
#' @param output If TRUE produces summary of model input and model output
#' @param title Optional title for the MARK analysis output
#' @param invisible if TRUE, exectution of MARK.EXE is hidden from view
#' @param adjust if TRUE, adjusts number of parameters (npar) to number of
#' columns in design matrix, modifies AIC and records both
#' @param se if TRUE, se and confidence intervals are shown in summary sent to
#' screen
#' @param filename base filename for files created by MARK.EXE. Files are named
#' filename.*.
#' @param prefix base filename prefix for files created by MARK.EXE; the files
#' are named prefixnnn.*
#' @param default.fixed if TRUE, real parameters for which the design data have
#' been deleted are fixed to default values
#' @param silent if TRUE, errors that are encountered are suppressed
#' @param retry number of reanalyses to perform with new starting values when
#' one or more parameters are singular
#' @param realvcv if TRUE the vcv matrix of the real parameters is extracted
#' and stored in the model results
#' @param external if TRUE the mark object is saved externally rather than in
#' the workspace; the filename is kept in its place
#' @param threads number of cpus to use with mark.exe if positive or number of cpus to remain idle if negative
#' @return model: MARK model object with the base filename stored in
#' \code{output} and the extracted \code{results} from the output file appended
#' onto list; see \code{\link{mark}} for a detailed description of a
#' \code{mark} object.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{make.mark.model}}, \code{\link{run.models}},
#' \code{\link{extract.mark.output}}, \code{\link{adjust.parameter.count}},
#' \code{\link{mark}}, \code{\link{cleanup}}
#' @keywords model
#' @examples
#' \dontrun{
#' # The following example will not run because the data are not included in the
#' # examples.  It illustrates the use of rerun.mark with mark.wrapper.  With this
#' # particular data set the POPAN models were having difficulty converging.  After
#' # running the set of models using mark.wrapper and looking at the results it
#' # was clear that in several instances the model did not converge.  This is easiest
#' # to discern by comparing nested models in the model.table.  If one model 
#' # is nested within another,then the deviance of the model with more
#' # parameters should be as good or better than the smaller model.  If that 
#' # is not the case then the model that converged can be used for initial 
#' # values in a call to rerun.mark for the model that did not converge.
#' #
#' 
#' do.nat=function()
#' {
#' Phi.ageclass=list(formula=~ageclass)
#' Phi.dot=list(formula=~1)
#' p.area=list(formula=~area)
#' p.timebin.plus.area=list(formula=~timebin+area)
#' p.timebin.x.area=list(formula=~-1+timebin:area)
#' pent.ageclass=list(formula=~ageclass)
#' pent.ageclass.plus.EN=list(formula=~ageclass+EN)
#' pent.ageclass.plus.diffEN=list(formula=~ageclass+EN92+EN97+EN02)
#' cml=create.model.list("POPAN")
#' nat=mark.wrapper(cml,data=zc.proc,ddl=zc.ddl,
#'   invisible=FALSE,initial=1,retry=2)
#' return(nat)
#' }
#' nat=do.nat()
#' # model list
#' #            Phi                   p                      pent
#' #1  Phi.ageclass              p.area             pent.ageclass
#' #2  Phi.ageclass              p.area pent.ageclass.plus.diffEN
#' #3  Phi.ageclass              p.area     pent.ageclass.plus.EN
#' #4  Phi.ageclass p.timebin.plus.area             pent.ageclass
#' #5  Phi.ageclass p.timebin.plus.area pent.ageclass.plus.diffEN
#' #6  Phi.ageclass p.timebin.plus.area     pent.ageclass.plus.EN
#' #7  Phi.ageclass    p.timebin.x.area             pent.ageclass
#' #8  Phi.ageclass    p.timebin.x.area pent.ageclass.plus.diffEN
#' #9  Phi.ageclass    p.timebin.x.area     pent.ageclass.plus.EN
#' #10      Phi.dot              p.area             pent.ageclass
#' #11      Phi.dot              p.area pent.ageclass.plus.diffEN
#' #12      Phi.dot              p.area     pent.ageclass.plus.EN
#' #13      Phi.dot p.timebin.plus.area             pent.ageclass
#' #14      Phi.dot p.timebin.plus.area pent.ageclass.plus.diffEN
#' #15      Phi.dot p.timebin.plus.area     pent.ageclass.plus.EN
#' #16      Phi.dot    p.timebin.x.area             pent.ageclass
#' #17      Phi.dot    p.timebin.x.area pent.ageclass.plus.diffEN
#' #18      Phi.dot    p.timebin.x.area     pent.ageclass.plus.EN
#' #
#' # use model 9 as starting values for model 7
#' nat[[7]]= rerun.mark(nat[[7]],data=zc.proc,ddl=zc.ddl,initial=nat[[9]])
#' # use model 3 as starting values for model 1
#' nat[[1]]= rerun.mark(nat[[1]],data=zc.proc,ddl=zc.ddl,initial=nat[[3]])
#' # use model 14 as starting values for model 15
#' nat[[15]]= rerun.mark(nat[[15]],data=zc.proc,ddl=zc.ddl,initial=nat[[14]])
#' # use model 5 as starting values for model 6
#' nat[[6]]= rerun.mark(nat[[6]],data=zc.proc,ddl=zc.ddl,initial=nat[[5]])
#' # use model 10 as starting values for model 11
#' nat[[11]]= rerun.mark(nat[[11]],data=zc.proc,ddl=zc.ddl,initial=nat[[10]])
#' # use model 10 as starting values for model 12
#' nat[[12]]= rerun.mark(nat[[12]],data=zc.proc,ddl=zc.ddl,initial=nat[[10]])
#' # reconstruct model table with new results
#' nat$model.table=model.table(nat[1:18])
#' # show new model table
#' nat
#' }
rerun.mark <-
function(model,data,ddl,initial,output=TRUE,title="",invisible=TRUE,adjust=TRUE,se=FALSE,
 filename=NULL,prefix="mark",default.fixed=TRUE,silent=FALSE,retry=0,realvcv=FALSE,external=FALSE,threads=-1)
{
# -----------------------------------------------------------------------------------------------------------------------
# rerun.mark -  reruns previous mark model with different initial values
#
# Arguments:
#
#  model                - previously run mark model
#  data                 - processed dataframe
#  ddl                  - design data list which contains an element for each parameter type
#  initial              - vector of initial values for beta parameters
#  output               - if TRUE produces summary of model input and model output
#  invisible            - if TRUE, window for running MARK is hidden
#  adjust               - if TRUE, adjusts npar to # of cols in design matrix, modifies AIC and records both
#  se                   - if TRUE, se and confidence intervals are shown in summary sent to screen
#  filename             - base filename for MARK input and output files
#  prefix               - base filename prefix; default is "mark" for files named marknnn.*
#  default.fixed        - if TRUE, default fixed values are assigned to any parameters missing from the full design data
#  silent               - if TRUE, errors that are encountered are suppressed
#  retry                - number of reanalyses to perform with new starting values when one or more parameters are singular
#  realvcv              - if TRUE the vcv matrix of the real parameters is extracted and stored in the model results
#  external             - if TRUE the mark object is saved externally rather than in the workspace; the filename is kept in its place
#
#  Value: 
#
#  model - a MARK object model containing output and extracted results
#
#  Functions used: make.mark.model, run.mark.model, summary.mark
# 
# -----------------------------------------------------------------------------------------------------------------------
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
# 
simplify=TRUE
model=load.model(model)
if(is.null(data$data))
   stop("\nMust specify processed dataframe\n")
#
# If the design data have not been specified, stop
#
if(is.null(ddl))
   stop("\nMust specify design data list\n")
#
#  Assign model.parameters
#
model.parameters=model$model.parameters
#
# Run model as many as times requested if needed
#
i=0
converge=FALSE
while(i<=retry & !converge)
{
#
# Remake the model with new initial values
#

   model<-make.mark.model(data,title=title,parameters=model.parameters,
          ddl=ddl,initial=initial,call=match.call(),default.fixed=default.fixed,
          model.name=model$model.name)
   model$model.parameters=model.parameters
#
# Summarize model input if output=TRUE
#
   if(output & i==1)
   {
     message("\n")
     print(summary(model))
   }
#
# Run model
#
   runmodel<-try(run.mark.model(model,invisible=invisible,adjust=adjust,filename=filename,prefix=prefix,realvcv=realvcv,threads=threads,),silent=silent)
   if(class(runmodel)[1]=="try-error")
     stop("\n\n********Following model failed to run :",model$model.name,"**********\n\n")
   else
   {
#
#  Check if any parameters are singular and if retry=TRUE, the refit model with
#  new initial values
#
      if(retry>0 && !is.null(runmodel$results$singular))
      {
         message("\nRe-running analysis with new starting values\n")
         i=i+1
         converge=FALSE
         initial=runmodel$results$beta$estimate
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
      message("\n")
      print(summary(runmodel,se=se))
   }
#
# Return fitted MARK model object 
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
