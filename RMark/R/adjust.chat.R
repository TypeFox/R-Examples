#' Adjust over-dispersion scale or a result value such as effective sample size
#' 
#' Adjust value of over-dispersion constant or another result value for a
#' collection of models which modifies model selection criterion and estimated
#' standard errors.
#' 
#' The value of \code{chat} is stored with the model object except when there
#' is no over-dispersion (\code{chat}=1). This function assigns a new value of
#' \code{chat} for the collection of models specified by \code{model.list}
#' and/or \code{type}.  The value of \code{chat} is used by
#' \code{\link{model.table}} for model selection in computing QAICc unless
#' \code{chat=1}.  It is also used in \code{\link{summary.mark}},
#' \code{\link{get.real}} and \code{\link{compute.real}} to adjust standard
#' errors and confidence intervals. Note that the standard errors and
#' confidence intervals in \code{results$beta},\code{results$beta.vcv}
#' \code{results$real}, \code{results$derived} and \code{results$derived.vcv}
#' are not modified and always assume \code{chat=1}.
#' 
#' It can also be used to modify a field in \code{model$results} such as
#' \code{n} which is ESS (effective sample size) from MARK output that is used
#' in AICc/QAICc calculations.
#' 
#' @aliases adjust.chat adjust.value
#' @usage adjust.value(field="n",value,model.list)
#'        adjust.chat(chat=1,model.list)
#' @param field Character string containing name of the field; either
#' \code{chat} or a field in \code{model$results} such as \code{n} for sample
#' size used in AICc or QAICc
#' @param value new value for field
#' @param chat Over-dispersion scale
#' @param model.list marklist created by the function
#' \code{\link{collect.models}} which has each model object and a
#' \code{model.table} at the end. For the entire collection of models each
#' \code{chat} is adjusted. If the argument type is specified the collected
#' models are limited to mark analyses with that specific type of model ("CJS")
#' @return model.list with all models given the new chat value and model.table
#' adjusted for chat values
#' @note See note in \code{\link{collect.models}}
#' @author Jeff Laake
#' @export adjust.chat adjust.value
#' @seealso \code{\link{model.table}}, \code{\link{summary.mark}},
#' \code{\link{get.real}} ,\code{\link{compute.real}}
#' @keywords utility
#' @examples
#' 
#' #
#' # The following are examples only to demonstrate selecting different 
#' # model sets for adjusting chat and showing model selection table. 
#' # It is not a realistic analysis.
#' #
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' do_example=function()
#' {
#' mod1=mark(dipper)
#' mod2=mark(dipper,model.parameters=list(Phi=list(formula=~time)))
#' mod3=mark(dipper,model="POPAN",initial=1)
#' cjs.results=collect.models(type="CJS")
#' cjs.results  # show model selection results for "CJS" models
#' }
#' cjs.results=do_example()
#' cjs.results
#' # adjust chat for all models to 2
#' cjs.results=adjust.chat(2,cjs.results) 
#' cjs.results
#' }
adjust.value <- function(field="n",value,model.list)
# ----------------------------------------------------------------------------------------
#
# adjust.value  - adjusts values of field in model list for a collection of models
#
# Arguments:
#
# field         - name of field chat, n etc
# value         - new value of field
# model.list    - a marklist created by collect.models
# 
# Value:  
#
#  model.list with all models given the new chat value and model.table adjusted for chat values
#
# Functions used: collect.model.names
#
# ----------------------------------------------------------------------------------------
{
#
# If no model list specified, collect models from parent.frame; if
# model.list is a list created by collect.models
#
if(!missing(model.list))
{
   if(class(model.list)=="marklist")
   {
      if(names(model.list)[length(model.list)]=="model.table")
         model.list=model.list[1:(length(model.list)-1)]
   }
}
else
{
   stop("A model.list must be given")
}
#
# For each model in the list store the new value of chat in it.
#
for (i in 1:length(model.list))
{
    if(is.character(model.list[[i]]))
    {
       model=load.model(model.list[[i]])
       if(field=="chat")
          model[field]=value
       else
          model[["results"]][field]=value       
       save(model,file=model.list[[i]])
    }
    else
       if(field=="chat")
          model.list[[i]][field]=value
       else
          model.list[[i]][["results"]][field]=value       
}
#
# Next recreate model.table
#
model.list$model.table=model.table(model.list)
class(model.list)="marklist"
return(model.list)
}

adjust.chat <- function(chat=1,model.list)
  return(adjust.value("chat",chat,model.list))
  

