#' Create table of MARK model selection results
#' 
#' Constructs a table of model selection results for MARK analyses. The table
#' includes the formulas, model name, number of parameters, deviance, AICc,
#' DeltaAICc, model weight and residual deviance.  If chat>1 QAICc, QDeltaAICc
#' and QDeviance are used instead.
#' 
#' This function is used by \code{\link{collect.models}} to construct a table
#' of model selection results with the models that it collects; however it can
#' be called directly to construct the table.
#' 
#' @param model.list a vector of model names or a list created by the function
#' \code{\link{collect.models}} which has each model object and at the end a
#' \code{model.table} ; If nothing is specified then any mark object in the
#' workspace is collected for the table.  If \code{type} is specified all
#' analyses in parent frame(\code{pf}) of that \code{type} of model are used.
#' If specified set of models are of conflicting types or of different data
#' sets then an error is issued unless ignore=TRUE
#' @param type type of model (eg "CJS")
#' @param sort if true sorts models by criterion
#' @param adjust if TRUE adjusts # of parameters to # of cols in design matrix
#' @param ignore if TRUE collects all models and ignores that they are from
#' different models
#' @param pf parent frame value; default=1 so it looks in calling frame of
#' model.table; used in other functions with pf=2 when functions are nested
#' two-deep
#' @param use.lnl display -2lnl instead of deviance
#' @param use.AIC use AIC instead of AICc
#' @param model.name if TRUE uses the model.name in each mark object which uses
#' formula notation.  If FALSE it uses the R names for the model obtained from
#' collect.model.names or names assigned to marklist elements
#' @return result.table - dataframe containing summary of models
#' \item{model.name}{name of fitted model} \item{parameter.name - an entry for
#' each parameter}{formula for parameter} \item{npar}{number of estimated
#' parameters} \item{AICc or QAICc}{AICc value or QAICc if chat>1}
#' \item{DeltaAICc or DeltaQAICc}{difference between AICc or QAICc value from
#' model with smallest value} \item{weight}{model weight based on
#' exp(-.5*DeltaAICc) or exp(-.5*QDeltaAICc)} \item{Deviance or
#' QDeviance}{residual deviance from saturated model}
#' \item{chat}{overdispersion constant if not 1}
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{collect.model.names}}, \code{\link{collect.models}}
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
#' run.dipper=function()
#' {
#' #
#' # Process data
#' #
#' dipper.processed=process.data(dipper,groups=("sex"))
#' #
#' # Create default design data
#' #
#' dipper.ddl=make.design.data(dipper.processed)
#' #
#' # Add Flood covariates for Phi and p that have different values
#' #
#' dipper.ddl$Phi$Flood=0
#' dipper.ddl$Phi$Flood[dipper.ddl$Phi$time==2 | dipper.ddl$Phi$time==3]=1
#' dipper.ddl$p$Flood=0
#' dipper.ddl$p$Flood[dipper.ddl$p$time==3]=1
#' #
#' #  Define range of models for Phi
#' #
#' Phi.dot=list(formula=~1)
#' Phi.time=list(formula=~time)
#' Phi.sex=list(formula=~sex)
#' Phi.sextime=list(formula=~sex+time)
#' Phi.sex.time=list(formula=~sex*time)
#' Phi.Flood=list(formula=~Flood)
#' #
#' #  Define range of models for p
#' #
#' p.dot=list(formula=~1)
#' p.time=list(formula=~time)
#' p.sex=list(formula=~sex)
#' p.sextime=list(formula=~sex+time)
#' p.sex.time=list(formula=~sex*time)
#' p.Flood=list(formula=~Flood)
#' #
#' # Return model table and list of models
#' #
#' cml=create.model.list("CJS")
#' return(mark.wrapper(cml,data=dipper.processed,ddl=dipper.ddl))
#' }
#' 
#' dipper.results=run.dipper()
#' dipper.results
#' dipper.results$model.table=model.table(dipper.results,model.name=FALSE)
#' dipper.results
#' #
#' # Compute matrices of model weights, number of parameters and Delta AICc values
#' #
#' model.weight.matrix=tapply(dipper.results$model.table$weight,
#'  list(dipper.results$model.table$Phi,dipper.results$model.table$p),mean)
#' model.npar.matrix=tapply(dipper.results$model.table$npar,
#'  list(dipper.results$model.table$Phi,dipper.results$model.table$p),mean)
#' model.DeltaAICc.matrix=tapply(dipper.results$model.table$DeltaAICc,
#'  list(dipper.results$model.table$p,dipper.results$model.table$Phi),mean)
#' #
#' # Output DeltaAICc as a tab-delimited text file that can be read into Excel 
#' # (to do that directly use RODBC or xlsreadwrite package for R)
#' #
#' write.table(model.DeltaAICc.matrix,"DipperDeltaAICc.txt",sep="\t")
#' }
model.table <-
function(model.list=NULL,type=NULL,sort=TRUE,adjust=TRUE,ignore=TRUE,pf=1,
              use.lnl=FALSE,use.AIC=FALSE,model.name=TRUE)
# ----------------------------------------------------------------------------------------
#
# model.table  - constructs a table of mark analyses for model selection. 
#                The table includes the model name, # of parameters, deviance, AICc
#                DeltaAICc and model weight.  If chat>1 QAICc is used instead
#
# Arguments:
#
# model.list    - a vector of model names or a list created by the function collect.models
#                 which has a model.table followed by each model object;
#                 If nothing is specified then any mark object in the workspace is collected for the
#                 table.  If type is specified all analyses in calling 
#                 frame of type model are used. If specified set of models are of conflicting 
#                 types or of different data sets then an error is issued unless ignore=TRUE
# type          - type of model (eg "CJS")
# sort          - if true sorts models by criterion
# adjust 	      - if TRUE uses adjusted # of parameters and adjusted AICc
#                   10 Jan 06 default changed to true
# ignore 	      - if TRUE collects all models and ignores that they are from different models
# pf            - parent frame value
# use.lnl       - display -2lnl instead of deviance
# use.AIC       - use AIC instead of AICc
# model.name    - if TRUE uses the model.name in each mark object which uses formula notation.
#                 If FALSE it uses the R names for the model obtained from collect.model.names
#                 or names assigned to marklist elements.
#
# Value:  
#
#  result.table - dataframe containing summary
#
# Functions used: collect.model.names
#
# ----------------------------------------------------------------------------------------
{
#
# If no model list specified, collect models from parent.frame
#
deviance=NULL
qdeviance=NULL
warned=FALSE
if(is.list(model.list))
{
   model.table=model.list
   model.list=names(model.list)
   if(model.list[length(model.list)]=="model.table")
      model.list=model.list[1:(length(model.list)-1)]
   amodeltable=TRUE
}
else
{
   amodeltable=FALSE
   lx = ls(envir = parent.frame(pf))
   if (is.null(model.list))
       model.list = collect.model.names(lx, type)
}
#
# For each model in the list, extract the relevant bits and store in the result.table
#
result.table=list()
model.numbers=NULL
#chat.values=rep(1,length(model.list))
chat.values=NULL
first=TRUE
for(i in 1:length(model.list))
{
    if (!amodeltable)
    {
       if(!(model.list[i] %in% lx))
         stop(paste(model.list[i], "model not found or model type not recognized\n"))
       model = eval(parse(text = model.list[i]), envir = parent.frame(pf))
       model=load.model(model)
    }
    else
      model = eval(parse(text = paste("load.model(",model.list[i],")")), envir = model.table)
   if(is.null(model$result))
   {
        message(paste("Model ",i,": ",model.list[i]," failed to run.\n",sep=""))
        next
   }
   param.names=names(setup.parameters(model$model))
   if(is.null(model$chat))
      chat=1
   else
      chat=model$chat
   chat.values=c(chat.values,chat)
#
#  If this is the first model, save the type and the effective sample size
#
   if(first)
   {
      savetype=model$model
      saveESS=model$result$n
	  first=FALSE
   }
#                     
#  If the models don't agree on type or ESS then stop
#
   else
   {
      if(savetype!=model$model)
      {
          if(ignore)
          {
             if(!warned)
             {
               warned=TRUE
               warning("Model list contains models of differing types\n")
               ncol=dim(result.table)[2]
               result.table=result.table[,(ncol-3):ncol]
             }
          }
          else
	     stop("Model list contains models of differing types")
      }
      if(saveESS!=model$result$n)stop("Model list contains models from different sets of data")
   }
#                                           
#  Restore unadjusted values if any and adjust=FALSE
#
   if(!adjust)
      if(!is.null(model$results$npar.unadjusted))
      {
         model$results$npar=model$results$npar.unadjusted
         model$results$AICc=model$results$AICc.unadjusted
      }
#
#  Store deviance, AICc, QAICc, and QDeviance
#
    if(use.lnl)
       deviance=c(deviance,model$result$lnl)
    else
       deviance=c(deviance,model$result$deviance)
    K=model$results$npar
    if(use.AIC)
       qaicc= model$results$lnl/chat + 2*K
    else
       qaicc= model$results$lnl/chat + 2*K + 2*K*(K+1)/(model$results$n-K-1)
    model.numbers=c(model.numbers,i)
    formulae=sapply(model$parameters,function(x){return(paste(x$formula,collapse=""))})
    if(use.AIC)
      AIC.value=model$result$lnl+2*K
    else
      AIC.value=model$results$lnl + 2*K + 2*K*(K+1)/(model$results$n-K-1)
    if(model.name) mn=model$model.name else mn=model.list[i]
    if(warned)
       result.table=rbind(result.table,data.frame(model=mn,npar=model$result$npar,AICc=AIC.value,QAICc=qaicc))
    else
       result.table=rbind(result.table,data.frame(t(formulae),model=mn,npar=model$result$npar,AICc=AIC.value,QAICc=qaicc))
    if(use.lnl)
       qdeviance=c(qdeviance,model$result$lnl)
    else
       qdeviance=c(qdeviance,model$result$deviance/chat)
}
row.names(result.table)=model.numbers
if(any(chat.values!=chat.values[1]))warning("Different chat values in collection of models\n")
anychat=any(chat.values!=1)
# 
# Compute model weight
#
if(!anychat)
{
   result.table$DeltaAICc=result.table$AICc-min(result.table$AICc)
   result.table$weight=exp(-.5*result.table$DeltaAICc)
   result.table$QAICc=NULL
}
else
{
   result.table$DeltaQAICc=result.table$QAICc-min(result.table$QAICc)
   result.table$weight=exp(-.5*result.table$DeltaQAICc)
   result.table$AICc=NULL
}
result.table$weight=result.table$weight/sum(result.table$weight)
if(!anychat)
   result.table$Deviance=deviance
else
{
   result.table$QDeviance=qdeviance
   result.table$chat=chat.values
}
#
# Exclude very tiny weights so output is easier to read
#
result.table$weight[result.table$weight<1e-15]=0
result.table$weight=result.table$weight/sum(result.table$weight)
#
# Sort if requested, otherwise return result.table in the order they were specified in the list
#
if(sort)
{
   if(!anychat)
      result.table=result.table[order(result.table$DeltaAICc),]
   else
      result.table=table=result.table[order(result.table$DeltaQAICc),]
}
#
#  Change names to match choices
#
   if(use.lnl)
   {
      names(result.table)[names(result.table)=="Deviance"]="Neg2LnL"
      names(result.table)[names(result.table)=="QDeviance"]="Neg2LnL"
   }
   if(use.AIC)
   {
      if(!anychat)
      {
         names(result.table)[names(result.table)=="AICc"]="AIC"
         names(result.table)[names(result.table)=="DeltaAICc"]="DeltaAIC"
      }
      else
      {
         names(result.table)[names(result.table)=="QAICc"]="QAIC"
         names(result.table)[names(result.table)=="DeltaQAICc"]="DeltaQAIC"
      }
   }
return(result.table)
}
