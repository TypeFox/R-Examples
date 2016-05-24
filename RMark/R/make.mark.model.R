#' Create a MARK model for analysis
#' 
#' Creates a MARK model object that contains a MARK input file with PIMS and
#' design matrix specific to the data and model structure and formulae
#' specified for the model parameters.
#' 
#' This function is called by \code{\link{mark}} to create the model but it can
#' be called directly to create but not run the model.  All of the arguments
#' have default values except for the first 2 that specify the processed data
#' list (\code{data}) and the design data list (\code{ddl}). If only these 2
#' arguments are specified default models are used for the parameters.  For
#' example, following with the example from \code{\link{process.data}} and
#' \code{\link{make.design.data}}, the default model can be created with:
#' 
#' \code{mymodel=make.mark.model(proc.example.data,ddl)}
#' 
#' The call to \code{make.mark.model} creates a model object but does not do
#' the analysis.  The function returns a list that includes several fields
#' including a design matrix and the MARK input file that will be used with
#' \code{MARK.EXE} to run the analysis from function
#' \code{\link{run.mark.model}}. The following shows the names of the list
#' elements in mymodel:
#' 
#' \preformatted{ names(mymodel) [1] "data" "model" "title" "model.name"
#' "links" [6] "mixtures" "call" "parameters" "input" "number.of.groups" [11]
#' "group.labels" "nocc" "begin.time" "covariates" "fixed" [16] "design.matrix"
#' "pims" "design.data" "strata.labels" "mlogit.list" [21] "simplify" }
#' 
#' The list is defined to be a mark object which is done by assigning a class
#' vector to the list.  The classes for an R object can be viewed with the
#' class function as shown below:
#' 
#' \preformatted{ class(mymodel) [1] "mark" "CJS" } Each MARK model has 2 class
#' values.  The first identifies it as a mark object and the second identifies
#' the type of mark analysis, which is the default "CJS" (recaptures only) in
#' this case.  The use of the class feature has advantages in using generic
#' functions and identifying types of objects.  An object of class \code{mark}
#' is defined in more detail in function \code{\link{mark}}.
#' 
#' To fit non-trivial models it is necessary to understand the remaining
#' calling arguments of \code{make.mark.model} and R formula notation. The
#' model formulae are specified with the calling argument \code{parameters}.
#' It uses a similar list structure as the \code{parameters} argument in
#' \code{\link{make.design.data}}. It expects to get a list with elements named
#' to match the parameters in the particular analysis (e.g., Phi and p in CJS)
#' and each list element is a list, so it is a list of lists).  For each
#' parameter, the possible list elements are \code{formula, link, fixed,
#' component, component.name, remove.intercept}. In addition, for closed
#' capture models and robust design model, the element \code{share} is included
#' in the list for p (capture probabilities) and GammaDoublePrime
#' (respectively) to indicate whether the model is shared (share=TRUE) or
#' not-shared (the default) (share=FALSE) with c (recapture probabilities) and
#' GammaPrime respectively.
#' 
#' \code{formula} specifies the model for the parameter using R formula
#' notation. An R formula is denoted with a \code{~} followed by variables in
#' an equation format possibly using the \code{+} , \code{*}, and \code{:}
#' operators.  For example, \code{~sex+age} is an additive model with the main
#' effects of \code{sex} and \code{age}. Whereas, \code{~sex*age} includes the
#' main effects and the interaction and it is equivalent to the formula
#' specified as \code{~sex+age+sex:age} where \code{sex:age} is the interaction
#' term.  The model \code{~age+sex:age} is slightly different in that the main
#' effect for \code{sex} is dropped which means that intercept of the
#' \code{age} effect is common among the sexes but the age pattern could vary
#' between the sexes.  The model \code{~sex*Age} which is equivalent to
#' \code{~sex + Age + sex:Age} has only 4 parameters and specifies a linear
#' trend with age and both the intercept and slope vary by sex.  One additional
#' operator that can be useful is \code{I()} which allows computation of
#' variables on the fly.  For example, the addition of the Agesq variable in
#' the design data (as described above) can be avoided by using the notation
#' \code{~Age + I(Age^2)} which specifies use of a linear and quadratic effect
#' for age.  Note that specifying the model \code{~age + I(age^2)} would be
#' incorrect and would create an error because \code{age} is a factor variable
#' whereas \code{Age} is not.
#' 
#' As an example, consider developing a model in which Phi varies by age and p
#' follows a linear trend over time. This model could be specified and run as
#' follows:
#' 
#' \preformatted{ p.Time=list(formula=~Time) Phi.age=list(formula=~age)
#' Model.p.Time.Phi.age=make.mark.model(proc.example.data,ddl,
#' parameters=list(Phi=Phi.age,p=p.Time))
#' Model.p.Time.Phi.age=run.mark.model(Model.p.Time.Phi.age) }
#' 
#' The first 2 commands define the p and Phi models that are used in the
#' \code{parameter} list in the call to \code{make.mark.model}. This is a good
#' approach for defining models because it clearly documents the models, the
#' definitions can then be used in many calls to \code{make.mark.model} and it
#' will allow a variety of models to be developed quickly by creating different
#' combinations of the parameter models.  Using the notation above with the
#' period separating the parameter name and the description (eg., p.Time) gives
#' the further advantage that all possible models can be developed quickly with
#' the functions \code{\link{create.model.list}} and
#' \code{\link{mark.wrapper}}.
#' 
#' Model formula can use any field in the design data and any individual
#' covariates defined in \code{data}.  The restrictions on individual
#' covariates that was present in versions before 1.3 have now been removed.
#' You can now use interactions of individual covariates with all design data
#' covariates and products of individual covariates.  You can specify
#' interactions of individual covariates and factor variables in the design
#' data with the formula notation.  For example, \code{~region*x1} describes a
#' model in which the slope of \code{x1} varies by \code{region}. Also,
#' \code{~time*x1} describes a model in which the slope for \code{x1} varied by
#' time; however, there would be only one value of the covariate per animal so
#' this is not a time varying covariate model.  Models with time varying
#' covariates are now more easily described with the improvements in version
#' 1.3 but there are restrictions on how the time varying individual covariates
#' are named.  An example would be a trap dependence model in which capture
#' probability on occasion i+1 depends on whether they were captured on
#' occasion i.  If there are n occasions in a CJS model, the 0/1 (not
#' caught/caught) for occasions 1 to n-1 would be n-1 individual covariates to
#' describe recapture probability for occasions 2 to n. For times 2 to n, a
#' design data field could be defined such that the variable timex is 1 if
#' time==x and 0 otherwise.  The time varying covariates must be named with a
#' time suffix on the base name of the covariate. In this example they would be
#' named as \code{x2,. . .,xn} and the model could be specified as \code{~time
#' + x} for time variation and a constant trap effect or as \code{~time +
#' time:x} for a trap effect that varied by time.  If in the
#' \code{\link{process.data}} call, the argument \code{begin.time} was set to
#' the year 1990, then the variables would have to be named x1991,x1992,...
#' because the first recapture occasion would be 1991.  Note that the times are
#' different for different parameters.  For example, survival is labeled based
#' on the beginning of the time interval which would be 1990 so the variables
#' should be named appropriately for the parameter model in which they will be
#' used.
#' 
#' In previous versions to handle time varying covariates with a constant
#' effect, it was necessary to use the component feature of the parameter
#' specification to be able to append one or more arbitrary columns to the
#' design matrix.  That is no longer required for individual covariates and the
#' component feature was removed in v2.0.8.
#' 
#' There are three other elements of the parameter list that can be useful on
#' occasion.  The first is \code{link} which specifies the link function for
#' transforming between the beta and real parameters.  The possible values are
#' "logit", "log", "identity" and "mlogit(*)" where * is a numeric identifier.
#' The "sin" link is not allowed because all models are specified using a
#' design matrix.  The typical default values are assigned to each parameter
#' (eg "logit" for probabilities, "log" for N, and "mlogit" for pent in POPAN),
#' so in most cases it will not be necessary to specify a link function.
#' 
#' The second is \code{fixed} which allows real parameters to be set at fixed
#' values.  The values for \code{fixed} can be either a single value or a list
#' with 5 alternate forms for ease in specifying the fixed parameters.
#' Specifying \code{fixed=value} will set all parameters of that type to the
#' specified value.  For example, \code{Phi=list(fixed=1)} will set all Phi to
#' 1.  This can be useful for situations like specifying F in the
#' Burnham/Barker models to all have the same value of 1.  Fixed values can
#' also be specified as a list in which values are specified for certain
#' indices, times, ages, cohorts, and groups.  The first 3 will be the most
#' useful.  The first list format is the most general and flexible but it
#' requires an understanding of the PIM structure and index numbers for the
#' real parameters.  For example,
#' 
#' \code{Phi=list(formula=~time, fixed=list(index=c(1,4,7),value=1))}
#' 
#' specifies Phi varying by time, but the real parameters 1,4,7 are set to 1.
#' The \code{value} field is either a single constant or its length must match
#' the number of indices.  For example,
#' 
#' \code{Phi=list(formula=~time, fixed=list(index=c(1,4,7),value=c(1,0,1)))}
#' 
#' sets real parameters 1 and 7 to 1 and real parameter 4 to 0.  Technically,
#' the index/value format for fixed is not wedded to the parameter type (i.e.,
#' values for p can be assigned within Phi list), but for the sake of clarity
#' they should be restricted to fixing real parameters associated with the
#' particular parameter type.  The \code{time} and \code{age} formats for fixed
#' will probably be the most useful.  The format fixed=list(time=x, value=y)
#' will set all real parameters (of that type) for time x to value y.  For
#' example,
#' 
#' \code{p=list(formula=~time,fixed=list(time=1986,value=1))}
#' 
#' sets up time varying capture probability but all values of p for 1986 are
#' set to 1.  This can be quite useful to set all values to 0 in years with no
#' sampling (e.g., \preformatted{fixed=list(time=c(1982,1984,1986), value=0)}).
#' The \code{age}, \code{cohort} and \code{group} formats work in a similar
#' fashion.  It is important to recognize that the value you specify for
#' \code{time}, \code{age}, \code{cohort} and \code{group} must match the
#' values in the design data list.  This is another reason to add binned fields
#' for age, time etc with \code{\link{add.design.data}} after creating the
#' default design data with \code{\link{make.design.data}}.  Also note that the
#' values for \code{time} and \code{cohort} are affected by the
#' \code{begin.time} argument specified in \code{\link{process.data}}.  Had I
#' not specified \code{begin.time=1980}, to set p in the last occasion (1986),
#' the specification would be
#' 
#' \code{p=list(formula=~time,fixed=list(time=7,value=1))}
#' 
#' because begin.time defaults to 1.  The advantage of the time-, age-, and
#' cohort- formats over the index-format is that it will work regardless of the
#' group definition which can easily be changed by changing the \code{groups}
#' argument in \code{\link{process.data}}.  The index-format will be dependent
#' on the group structure because the indexing of the PIMS will most likely
#' change with changes in the group structure.
#' 
#' Parameters can also be fixed at default values by deleting the specific rows
#' of the design data. See \code{\link{make.design.data}} and material below.
#' The default value for fixing parameters for deleted design data can be
#' changed with the \code{default=value} in the parameter list.
#' 
#' The final useful element of the parameter list is the
#' \code{remove.intercept} argument.  It is set to TRUE to forcefully remove
#' the intercept.  In R notation this can be done by specifiying the formula
#' notation ~-1+... but in formula with nested interactions of factor variables
#' and additive factor variables the -1 notation will not remove the intercept.
#' It will simply adjust the column definitions but will keep the same number
#' of columns and the model will be overparameterized.  The problem occurs with
#' nested factor variables like tostratum within stratum for multistrata
#' designs (see \code{\link{mstrata}}).  As shown in that example, you can
#' build a formula -1+stratum:tostratum to have transitions that are
#' stratum-specific.  If however you also want to add a sex effect and you
#' specify -1+sex+stratum:tostratum it will add 2 columns for sex labelled M
#' and F when in fact you only want to add one column because the intercept is
#' already contained within the stratum:tostratum term.  The argument
#' remove.intercept will forcefully remove the intercept but it needs to be
#' able to find a column with all 1's.  For example,
#' Psi=list(formula=~sex+stratum:tostratum,remove.intercept=TRUE) will work but
#' Psi=list(formula=~-1+sex+stratum:tostratum,remove.intercept=TRUE) will not
#' work.  Also, the -1 notation should be used when there is not an added
#' factor variable because
#' \preformatted{Psi=list(formula=~stratum:tostratum,remove.intercept=TRUE)}
#' will not work because while the stratum:tostratum effectively includes an
#' intercept it is equivalent to using an identity matrix and is not specified
#' as treatment contrast with one of the columns as all 1's.
#' 
#' The argument simplify determines whether the pims are simplified such that
#' only indices for unique and fixed real parameters are used.  For example,
#' with an all different PIM structure with CJS with K occasions there are
#' K*(K-1) real parameters for Phi and p.  However, if you use simplify=TRUE
#' with the default model of Phi(.)p(.), the pims are re-indexed to be 1 for
#' all the Phi's and 2 for all the p's because there are only 2 unique real
#' parameters for that model.  Using simplify can speed analysis markedly and
#' probably should always be used.  This was left as an argument only to test
#' that the simplification was working and produced the same likelihood and
#' real parameter estimates with and without simplification. It only adjust the
#' rows of the design matrix and not the columns.  There are some restrictions
#' for simplification.  Real parameters that are given a fixed value are
#' maintained in the design matrix although it does simplify amongst the fixed
#' parameters.  For example, if there are 50 real parameters all fixed to a
#' value of 1 and 30 all fixed to a value of 0, they are reduced to 2 real
#' parameters fixed to 1 and 0. Also, real parameters such as Psi in
#' Multistrata and pent in POPAN that use multinomial logits are not simplified
#' because they must maintain the structure created by the multinomial logit
#' link.  All other parameters in those models are simplified.  The only
#' downside of simplification is that the labels for real parameters in the
#' MARK output are unreliable because there is no longer a single label for the
#' real parameter because it may represent many different real parameters in
#' the all-different PIM structure.  This is not a problem with the labels in R
#' because the real parameter estimates are translated back to the
#' all-different PIM structure with the proper labels.
#' 
#' The argument \code{default.fixed} is related to deletion of design data (see
#' \code{\link{make.design.data}}).  If design data are deleted and
#' \code{default.fixed=T} the missing real parameters are fixed at a reasonable
#' default to represent structural "zeros".  For example, p is set to 0, Phi is
#' set to 1, pent is set to 0, etc.  For some parameters there are no
#' reasonable values (e.g., N in POPAN), so not all parameters will have
#' defaults.  If a parameter does not have a default or if
#' \code{default.fixed=F} then the row in the design matrix for that parameter
#' is all zeroes and its real value will depend on the link function.  For
#' example, with "logit" link the real parameter value will be 0.5 and for the
#' log link it will be 1.  As long as the inverse link is defined for 0 it will
#' not matter in those cases in which the real parameter is not used because it
#' represents data that do not exist.  For example, in a "CJS" model if initial
#' captures (releases) only occur every other occasion, but recaptures
#' (resightings) occurred every occasion, then every other cohort (row) in the
#' PIM would have no data.  Those rows (cohorts) could be deleted from the
#' design data and it would not matter if the real parameter was fixed.
#' However, for the case of a Jolly-Seber type model (eg POPAN or Pradel
#' models) in which the likelihood includes a probability for the leading
#' zeroes and first 1 in a capture history (a likelihood component for the
#' first capture of unmaked animals), and groups represent cohorts that enter
#' through time, you must fix the real parameters for the unused portion of the
#' PIM (ie for occasions prior to time of birth for the cohort) such that the
#' estimated probability of observing the structural 0 is 1.  This is easily
#' done by setting the pent (probability of entry) to 0 or by setting the
#' probability of capture to 0 or both.  In that case if
#' \code{default.fixed=F}, the probabilities for all those parameters would be
#' incorrectly set to 0.5 for p and something non-zero but not predetermined
#' for pent because of the multinomial logit.  Now it may be possible that the
#' model would correctly estimate these as 0 if the real parameters were kept
#' in the design, but we know what those values are in that case and they need
#' not be estimated. If it is acceptable to set \code{default.fixed=F}, the
#' functions such as \code{\link{summary.mark}} recognize the non-estimated
#' real parameters and they are not shown in the summaries because in essence
#' they do not exist.  If \code{default.fixed=T} the parameters are displayed
#' with their fixed value and for \code{summary.mark(mymodel,se=TRUE)}, they
#' are listed as "Fixed".
#' 
#' Missing design data does have implications for specifying formula but only
#' when interactions are specified.  With missing design data various factors
#' may not be fully crossed.  For example, with 2 factors A and B, each with 2
#' levels, the data are fully crossed if there are data with values A1&B1,
#' A1&B2, A2&B1 and A2&B2.  If data exist for each of the 4 combinations then
#' you can described the interaction model as ~A*B and it will estimate 4
#' values.  However, if data are missing for one of more of the 4 cells then
#' the "interaction" formula should be specified as ~-1+A:B or ~-1+B:A or
#' ~-1+A%in%B or ~-1+B%in%A to estimate all of the parameters represented by
#' the combinations with data.  An example of this could be a marking program
#' with multiple sites which resighted at all occasions but which only marked
#' at sites on alternating occasions.  In that case time is nested within site
#' and time-site interaction models would be specified as ~-1+time:site.
#' 
#' The argument \code{title} supplies a character string that is used to label
#' the output. The argument \code{model.name} is a descriptor for the model
#' used in the analysis.  The code constructs a model name from the formula
#' specified in the call (e.g., \code{Phi(~1)p(~time)}) but on occasion the
#' name may be too long or verbose, so it can be over-ridden with the
#' \code{model.name} argument.
#' 
#' The final argument \code{initial} can be used to provide initial estimates
#' for the beta parameters. It is either 1) a single starting value for each
#' parameter, 2) an unnamed vector of values (one for each parameter), 3) named
#' vector of values, or 4) the name of \code{mark} object that has already been
#' run. For cases 3 and 4, the code only uses appropriate initial beta
#' estimates in which the column names of the design matrix (for model) or
#' vector names match the column names in the design matrix of the model to be
#' run.  Any remaining beta parameters without an initial value specified are
#' assigned a 0 initial value.  If case 4 is used the models must have the same
#' number of rows in the design matrix and thus presumably the same structure.
#' As long as the vector elements are named (#3), the length of the
#' \code{initial} vector no longer needs to match the number of parameters in
#' the new model as long as the elements are named. The names can be retrieved
#' either from the column names of the design matrix or from
#' \code{rownames(x$results$beta)} where \code{x} is the name of the
#' \code{mark} object.
#' 
#' \code{options} is a character string that is tacked onto the \code{Proc
#' Estimate} statement for the MARK .inp file.  It can be used to request
#' options such as NoStandDM (to not standardize the design matrix) or
#' SIMANNEAL (to request use of the simulated annealing optimization method) or
#' any existing or new options that can be set on the estimate proc.
#' 
#' @param data Data list resulting from function \code{\link{process.data}}
#' @param ddl Design data list from function \code{\link{make.design.data}}
#' @param parameters List of parameter formula specifications
#' @param title Title for the analysis (optional)
#' @param model.name Model name to override default name (optional)
#' @param initial Vector of named or unnamed initial values for beta parameters
#' or previously run model (optional)
#' @param call Pass function call when this function is called from another
#' function (e.g.\code{\link{mark}}) (internal use)
#' @param default.fixed if TRUE, real parameters for which the design data have
#' been deleted are fixed to default values
#' @param options character string of options for Proc Estimate statement in
#' MARK .inp file
#' @param profile.int if TRUE, requests MARK to compute profile intervals
#' @param chat pre-specified value for chat used by MARK in calculations of
#' model output
#' @param simplify if FALSE, does not simplify PIM structure
#' @param input.links specifies set of link functions for parameters with non-simplified structure
#' @param parm.specific if TRUE, forces a link to be specified for each parameter
#' @param mlogit0 if TRUE, any real parameter that is fixed to 0 and has an mlogit link will 
#' have its link changed to logit so it can be simplified
#' @param hessian if TRUE specifies to MARK to use hessian rather than second partial matrix
#' @param accumulate if TRUE accumulate like data values into frequencies
#' @param icvalues numeric vector of individual covariate values for computation of real values
#' @param wrap if TRUE, data lines are wrapped to be length 80; if length of a row is not a 
#'   problem set to FALSE and it will run faster
#' @return model: a MARK object except for the elements \code{output} and
#' \code{results}. See \code{\link{mark}} for a detailed description of the
#' list contents.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{process.data}},\code{\link{make.design.data}},
#' \code{\link{run.mark.model}} \code{\link{mark}}
#' @keywords model
#' @examples
#' 
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(dipper)
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
#' Phidot=list(formula=~1)
#' Phitime=list(formula=~time)
#' PhiFlood=list(formula=~Flood)
#' #
#' #  Define range of models for p
#' #
#' pdot=list(formula=~1)
#' ptime=list(formula=~time)
#' #
#' # Make assortment of models
#' #
#' dipper.phidot.pdot=make.mark.model(dipper.processed,dipper.ddl,
#'   parameters=list(Phi=Phidot,p=pdot))
#' dipper.phidot.ptime=make.mark.model(dipper.processed,dipper.ddl,
#'   parameters=list(Phi=Phidot,p=ptime))
#' dipper.phiFlood.pdot=make.mark.model(dipper.processed,dipper.ddl,
#'   parameters=list(Phi=PhiFlood, p=pdot))
#' }
make.mark.model <-
function(data,ddl,parameters=list(),title="",model.name=NULL,initial=NULL,call=NULL,
		default.fixed=TRUE,options=NULL,profile.int=FALSE,chat=NULL,simplify=TRUE,
		input.links=NULL,parm.specific=FALSE,mlogit0=FALSE,hessian=FALSE,accumulate=TRUE,
		icvalues=NULL,wrap=TRUE)
{

#  *******************  INTERNAL FUNCTIONS    *********************************
#
#  print.pim: prints pim file to outfile for use in constructing MARK input
#
"print.pim" <-
function(pim,outfile)
{
#
# Arguments:
#
# pim     - pim matrix
# outfile - name of output file to write pim
#
# Value: None
#
# Define internal function "xp" that pastes non-zero values together with
# intervening spaces
#
xp=function(x){paste(x[x>0],collapse=" ")}
#
# For each row in the pim, apply xp to create a vector of concatenated values
# and then paste a ";" to end of each value
#
if(is.matrix(pim))
   strings=paste(apply(pim,1,xp),";")
else
   strings=paste(paste(pim,collapse=" "),";")
#
# Output strings
#
write(strings,outfile,append=TRUE)
return(NULL)
}
#
# spell: changes capitalization on links so it will be acceptable to mark interface
#
"spell" <- function(links)
{
  newlinks=links
  newlinks[newlinks=="logit"]="Logit"
  newlinks[newlinks=="mlogit"]="MLogit"
  newlinks[newlinks=="log"]="Log"
  newlinks[newlinks=="loglog"]="LogLog"
  newlinks[newlinks=="cloglog"]="CLogLog"
  newlinks[newlinks=="identity"]="Identity"
  newlinks[newlinks=="sin"]="Sin"
  return(newlinks)
}
#
#  realign.pims: realigns pim values to represent structure of design matrix
#
"realign.pims" <-
function(model){
#
#  Arguments:
#
#  model - a mark model object that has been created by make.mark.model
#
#  Value:
#
#  new.indices - a vector of new indices for the old PIM values.  The old
#                PIM values are 1:length(new.indices) and the new index is
#                the corresponding value.  For example, new.indices=c(1,1,2,2)
#                renumbers the structure 1,2,3,4 such that 1,2 are now 1
#                and 3,4 are now 2.
#
#  Get all the unique rows in the design matrix and paste all the values
#  in each row together.
#
   uniquevals=apply(unique(cbind(model$design.matrix,model$links)),1,paste,collapse="")
#
#  Get all the rows in the design matrix and paste all the values
#  in each row together.
#
   allvals=apply(cbind(model$design.matrix,model$links),1,paste,collapse="")
#
#  Find the corresponding sets of indices by matching allvals into uniquevals
#
   new.indices=match(allvals, uniquevals)
#
#  Next cope with fixed real parameters; first determine the unique fixed values
#
   uniquefixedvalues=unique(model$fixed$value)
#
#  Now create a parameter index for each of the unique fixed real parameters
#  assuming that they are all different from parameters defined by design
#  matrix (ie add onto max of uniqueindices),
#
   uniqueindices=match(model$fixed$value,uniquefixedvalues)+max(new.indices)
#
#  Assign these new indices to their position in the original PIM set
#
   new.indices[model$fixed$index]=uniqueindices
#
#  Some may overlap and others are new, so they need to be renumbered once again
#  eliminating extra ones
#
   new.indices=match(new.indices,sort(unique(new.indices)))
#
#  Simplification cannot occur for parameters with an mlogit link, so these
#  must be given new indices
#     Following lines force non-simplification for mlogit parameters
   mlogit.parameters=substr(model$links,1,6)=="mlogit" | substr(model$links,1,6)=="MLogit"
   new.indices[mlogit.parameters]=(1:sum(as.numeric(mlogit.parameters)))+max(new.indices)
   new.indices=match(new.indices,sort(unique(new.indices)))
#
#   This code that was added in v1.4.5 hsa been removed for v1.6.6 once it became
#   clear that you cannot simplify mlogit parameters even in the mlogit set
#   because the summation is only restricted based on the unique real parameters.
#   Some addtional code in make.mark.model that creates the mlogit.structure
#   can also probably be eliminated .
#
#   mlogit.parameters=substr(model$links,1,6)=="mlogit"
#   if(any(mlogit.parameters))
#   {
#      new.indices[mlogit.parameters]=max(new.indices[mlogit.parameters])*(mlogit.list$structure-1)+
#                                     new.indices[mlogit.parameters] + max(new.indices)
#
#    The following code was added to adjust for different model structures within
#    a stratum across tostratum.  Without it the mlogit structure doesn't not work
#    for models with sum tostratum interactions.
#
#      for (i in 1:max(mlogit.list$structure))
#      {
#         mm=matrix(new.indices[mlogit.parameters][mlogit.list$structure==i],ncol=mlogit.list$ncol)
#         umm = unique(mm)
#         rownames(umm) = 1:nrow(umm)
#         if(nrow(umm)>1)
#         {
#            unvals=apply(umm,1,paste,collapse="")
#            mmvals=apply(mm,1,paste,collapse="")
#            umm=max(new.indices)+(col(umm)-1)*nrow(umm)+row(umm)
#            mm=umm[as.numeric(names(unvals)[match(mmvals,unvals)]),]
#            new.indices[mlogit.parameters][mlogit.list$structure==i]=as.vector(mm)
#         }
#     }
#      new.indices=match(new.indices,sort(unique(new.indices)))
#  }
   return(new.indices)
}
#
#  renumber.pims: using the vector of new indices that match the old
#                 structure to the new structure, change the values in the
#                 pim argument.  The way this is done depends on whether it is
#                 a square or triangular PIM.
#

"renumber.pims" <-
function(pim,newlist,type){
if(type%in%c("Triang","STriang"))
{
#   pim=t(pim)
	pim[pim!=0]=newlist[pim]
#    pim[lower.tri(pim,TRUE)]=newlist[pim]
#   return(t(pim))
    return(pim)
} else
   return(newlist[pim])
}
#check.mlogits=function(model)
#{
#   mlogit.parameters=substr(model$links,1,6)=="mlogit" | substr(model$links,1,6)=="MLogit"
#   if(any(mlogit.parameters))
#   {
#      for (i in 1:max(mlogit.list$structure))
#      {
#         mm=matrix(model$links[mlogit.parameters][mlogit.list$structure==i],ncol=mlogit.list$ncol)
#         if(!all(apply(mm,1,function(x)return(length(unique(x))==1))))
#         {
#            cat("\n Error in mlogit assignment. This should not happen. Contact maintainer.\n")
#            print(mm)
#            stop()
#         }
#      }
#   }
#   invisible()
#}

"pim.header"<- function(group,param.name,parameters,ncol,stratum,tostratum,strata.labels,mixtures,session=NULL,socc=NULL,bracket=FALSE)
{
  if(!is.null(stratum)&length(strata.labels)>0)
  {
	  if(bracket)stratum.designation=""	 
	  if(!is.null(tostratum))
	  {
		  if(bracket)
			  param.name=paste(param.name,"[",stratum,",",tostratum,"]",sep="")
		  else
			  stratum.designation=paste(stratum,"to",tostratum)
#		  stratum.designation=paste(strata.labels[stratum],"to",strata.labels[tostratum])
	  }
	  else
	  {
		  if(bracket)
			  param.name=paste(param.name,"[",stratum,"]",sep="")
		  else
			  stratum.designation= paste(stratum,":Stratum",stratum,sep="")
#		  stratum.designation= paste(strata.labels[stratum],":Stratum",stratum,sep="")
	  }
  }
  else
     stratum.designation=""
  if(is.null(session))
     session.designation=""
  else
     if(is.null(socc))
	    session.designation=paste("Session",session)
	 else
		 if(!socc)
			 session.designation=paste("Primary",session)
		 else
		     session.designation=paste("Sampling Occasion",session)
 if (parameters$type == "Triang")
         string = paste(paste("group=", group,sep=""), param.name, stratum.designation, session.designation, 
                           paste(" rows=",ncol," cols=",ncol,sep=""), parameters$type, ";")
 else
	 if (parameters$type == "STriang")
				   string = paste(paste("group=", group,sep=""), param.name, stratum.designation, session.designation, 
						   paste(" rows=",(mixtures+parameters$rows)*ncol," cols=",ncol,sep=""), parameters$type, ";")
	else
      if(mixtures==1)
          string=paste(paste("group=",group,sep=""),param.name,stratum.designation, session.designation, 
                           paste(" rows=1"," cols=",ncol,sep=""),parameters$type,";")
      else
          if(!is.null(parameters$mix) && parameters$mix)
            string=paste(paste("group=",group,sep=""),param.name,stratum.designation,session.designation, 
              paste(" rows=",mixtures+parameters$rows," cols=",ncol,sep=""),parameters$type,";")
          else
            string=paste(paste("group=",group,sep=""),param.name,stratum.designation,session.designation, 
                          paste(" rows=1"," cols=",ncol,sep=""),parameters$type,";")
return(string)
}
"simplify.pim.structure" <-
function(model)
{
#
# simplify.pim.structure: renumbers PIMS to represent model structure that
#                         was created with the formula. It takes the input
#                         for MARK created in the model by make.mark.model
#                         with the formulas and simplies the PIM structure
#                         represented by the unique rows in the design matrix.
#                         It recreates the new input for MARK to reflect the
#                         change and adds it and some other fields to model
#                         for the pim translation.
#
# Arguments:
#
#  model - a mark model object that has been created by make.mark.model
#
#
# Value:
#
#  model - same mark model object with added list elements simplify
#          and a rewritten input object for MARK
#
#
#
#  Beginning of simplify.pim.structure function; it recreates input for
#  MARK and uses an outfile like make.mark.model
#
tempfilename=tempfile("markxxx",tmpdir=getwd(),fileext=".tmp")
outfile=file(tempfilename,open="wt")
#
# Use realign.pims to simplify PIM structure
#
new.indices=realign.pims(model)
#
# Copy first portion of the MARK input file because it will be unchanged by
# the simplification
#
input=model$input[1:grep("model=\\{",model$input)]
writeLines(input,outfile)
#
#  If there are fixed real parameters then these need to be included in the
#  design matrix if they are not already done so by the formula.
#
if(!is.null(model$fixed))
{
#
#  Get the unique new indices,values from the original indices for the fixed parameters
#  in model$fixed$index
#
   fixed.parms=unique(data.frame(index=new.indices[model$fixed$index],value=model$fixed$value))
   fixed=NULL
   num.fixed=dim(fixed.parms)[1]
#
#  For each fixed real parameter write out the input strings for MARK to assign
#  the number of fixed values and each fixed index (parm) and its value.
#
   for (i in 1:num.fixed)
      fixed = c(fixed, paste("parm(", fixed.parms$index[i], ")=", fixed.parms$value[i],
                    sep = ""))
   string = paste("fixed=", num.fixed, ";",sep="")
   write(string, file = outfile, append = TRUE)
   write(paste(fixed, ";"), file = outfile, append = TRUE)
}
#
#  Assign some values from the model that are used below
#
parameters=model$parameters
param.names=sub("DoublePrime","''",names(parameters))
param.names=sub("Prime","'",param.names)
number.of.groups=model$number.of.groups
nocc=model$nocc
if(is.null(model$mixtures))
   mixtures=1
else
   mixtures=model$mixtures
#
#  For each type of parameter, output the new PIM structure; This largely
#  follows code in make.mark.model except it uses renumber.pims and print.pim
#  to renumber and print out the pim structure.
#
if(model$data$model=="RDMSOccRepro")
	bracket=TRUE
else
	bracket=FALSE
for (i in 1:length(parameters)) {
     for (j in 1:length(model$pims[[i]]))
     {
         ncol = dim(model$pims[[i]][[j]]$pim)[2]
         string=pim.header(pim[[i]][[j]]$group,param.names[i],parameters[[i]],
                   ncol,model$pims[[i]][[j]]$stratum,model$pims[[i]][[j]]$tostratum,model$strata.labels,
				   mixtures,model$pims[[i]][[j]]$session,parameters[[i]]$socc,bracket=bracket)
         write(string, file = outfile, append = TRUE)
         if(parameters[[i]]$type %in% c("Triang","STriang"))
         {
            newpim=renumber.pims(model$pims[[names(model$parameters)[i]]][[j]]$pim,new.indices,parameters[[i]]$type)
            print.pim(newpim,outfile)
         }
         else
         {
             nmix=1
             if(mixtures>1)
               if(!is.null(parameters[[i]]$mix) && parameters[[i]]$mix)
                   nmix=mixtures+parameters[[i]]$rows
             for(k in 1:nmix)
             {
               newpim=renumber.pims(model$pims[[names(model$parameters)[i]]][[j]]$pim[k,],new.indices,parameters[[i]]$type)
               print.pim(newpim,outfile)
             }
         }
     }
}
#
# Next compute the new simplified design matrix.  rownums is the row numbers (indices)
# from the original design matrix but there is only one for each of the new
# parameters with indices 1:length(new.indices).  The new design matrix
# (complete.design.matrix) is obtained by subsetting the rows from the
# original design matrix matching rownums.  This is done using row.names to be
# able to use subset so it will always yield a dataframe.  Using indices for
# row numbers can result in a vector if there is only a single beta.
rownums=match(1:length(unique(new.indices)),new.indices)
# 22-Aug-05; change to use [rownums,] to accomodate fixed parameters
# 1 feb 06; modified to cope with single element selected
if(length(rownums)==1)
	complete.design.matrix=subset(model$design.matrix,1:dim(model$design.matrix)[1]%in%rownums)
else
	complete.design.matrix=model$design.matrix[rownums,,drop=FALSE]
# for any fixed parameter set row to all 0s
if(!is.null(model$fixed))
	for (i in 1:num.fixed)
	complete.design.matrix[fixed.parms$index[i],]="0"
					#
# Find any columns that are all 0; left from mlogit0 fix
#
dm=complete.design.matrix
allzero=apply(dm,2,function(x) all(x=="0"))
complete.design.matrix=dm[,!allzero,drop=FALSE]
#
# Look for setting of initial values in the input file; if found write them 
# exclude columns that are now all 0s.
#
if(length(grep("XXXinitialXXX ",model$input))!=0)
{
	initial=strsplit(model$input[grep("XXXinitialXXX ", model$input)]," ")[[1]]
	initial=initial[-c(1:2,length(initial))]
	if(any(allzero)) initial=initial[!allzero]
	string=paste("initial ",paste(initial,collapse=" ")," ;",collapse=" ")
    write(string, file = outfile, append = TRUE)
}
#
#  If profile intervals requested write out needed statements
#

if(is.null(model$chat))chat=1
if(is.numeric(model$profile.int))
{
   if(any(!model$profile.int%in%new.indices)) 
      stop(paste("Profile interval argument requests values not in beta: 1 to ",
                      ncol(model$design.matrix),"\n"))
      string=paste(paste("Profile Intervals chat=",format(chat,digits=5),sep=""),
                          paste(model$profile.int,collapse=" ")) 
      write(paste(string,";",sep=""),file=outfile,append=TRUE)   
}
else
{
    if(model$profile.int)
    {
      string=paste(paste("Profile Intervals chat=",format(chat,digits=5),sep=""),
                            paste(unique(new.indices),collapse=" "))
      write(paste(string,";",sep=""),file=outfile,append=TRUE)      
    } 
}
# 10 Jan 06; change to accomodate S(.) with known fate where design matrix can
# become a single element with simplification
#if(is.vector(complete.design.matrix))
#{
#   complete.design.matrix=as.matrix(complete.design.matrix)
#   row.names(complete.design.matrix)=row.names(model$design.matrix)[1]
#}
# if icvalues not null, write out values to input file
#
if(!is.null(icvalues))
{
	string = paste("icvalues=", paste(icvalues,collapse=","), ";",sep="")
	write(string, file = outfile, append = TRUE)		
}
#
# Write out the design matrix into the MARK input file
#
string = paste("design matrix constraints=", dim(complete.design.matrix)[1],
        " covariates=", dim(complete.design.matrix)[2], ";",sep="")
write(string, file = outfile, append = TRUE)
write.table(complete.design.matrix, file = outfile, eol = ";\n",
       sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
#
# If there is a link specification for models that use different links for
# each real parameter, write those out now shifting them for translation of
# indices.
#
# 11 Jan 06; modified code for new format of links input
#  1 feb 06; fixed links code to work with simplification - not done previously
#  6 Apr 06; added code to deal with simplification within mlogit parameters
if (!is.null(model$links))
{
   if(length(model$links)>1)
   {
     newlinks=model$links[match(unique(new.indices), new.indices)][order(unique(new.indices))]
     xi=grep("mlogit(",newlinks,fixed=TRUE)
     logit.numbers=as.numeric(gsub(")","",gsub("mlogit(","",newlinks[xi],fixed=TRUE),fixed=TRUE))
     logit.numbers=match(logit.numbers,sort(unique(logit.numbers)))
     newlinks[xi]=paste("mlogit(",logit.numbers,")",sep="")
     write(paste("links=",length(unique(new.indices)),";",sep=""),file=outfile,append=TRUE)
     write(paste(spell(newlinks),";",sep=""), file = outfile, append = TRUE)
   }
   else
     newlinks=NULL
}
else
   newlinks=NULL
#
# Write out labels for betas
#
if(!is.matrix(complete.design.matrix))
   complete.design.matrix=as.matrix(complete.design.matrix)
string = paste("blabel(", 1:dim(complete.design.matrix)[2], ")=", colnames(complete.design.matrix),";",sep="")
write(string, file = outfile, append = TRUE)
#
# Write out labels for real parameters
#
rnames = rep("", dim(complete.design.matrix)[1])
ipos = 0
string = paste("rlabel(", 1:dim(complete.design.matrix)[1], ")=", row.names(complete.design.matrix), ";",sep="")
write(string, file = outfile, append = TRUE)
#
#  Complete with stop statement; then read the outfile into the input vector to
#  store in the model object.  delete the output file and add the fields to the
#  model object and return it.
#
write("proc stop;", file = outfile, append = TRUE)
close(outfile)
outfile=file(tempfilename,open="rt")
text = readLines(outfile)
close(outfile)
unlink(tempfilename)
model$input=text
if(!is.null(newlinks))
   model$simplify=list(design.matrix=complete.design.matrix,pim.translation=new.indices,links=newlinks)
else
   model$simplify=list(design.matrix=complete.design.matrix,pim.translation=new.indices)
return(model)
}

#
# Create internal function to create a pim
#
create.pim=function(nocc,parameters,npar,mixtures)
{
    ncol=nocc+parameters$num
    mat=NULL
    if(parameters$type%in%c("Triang","STriang"))
    {
		nmix=1
		if(mixtures>1)
			if(!is.null(parameters$mix)&&parameters$mix)
				nmix=mixtures+parameters$rows
		for (j in 1:nmix)
		{
		   ncol=nocc+parameters$num
		   for(k in 1:(nocc+parameters$num))
           {
               if(parameters$pim.type=="all")
               {
                   mat=rbind(mat,c(rep(0,k-1),npar:(npar+ncol-1)))
                   npar=npar+ncol
               }
               else
               {
                  if(parameters$pim.type=="time")
                      mat=rbind(mat,c(rep(0,k-1),(npar+k-1):(npar+k-1+ncol-1)))
                  else
					  if(parameters$pim.type=="age")
						  mat=rbind(mat,c(rep(0,k-1),npar:(npar+ncol-1)))
					  else
						  mat=rbind(mat,c(rep(0,k-1),rep(npar,ncol)))
               }
               ncol=ncol-1
		   }
        }
#        if(parameters$pim.type!="all")
#            npar=max(mat)+1
   }
   else
   {
        nmix=1
        if(mixtures>1)
            if(!is.null(parameters$mix)&&parameters$mix)
                nmix=mixtures+parameters$rows
        for(k in 1:nmix)
        {
            mat=rbind(mat,npar:(npar+ncol-1))
            npar=npar+ncol
        }
   }
return(mat)
}
#
# Creates time-dependent covariates for age of nest if there
# is a field called AgeDay1 in the data.  Then the field NestAge
# can be used in the formula for S.
#
create.agenest.var=function(data,init.agevar,time.intervals)
{
   nocc=length(time.intervals)
   age.mat=matrix(data[,init.agevar],nrow=dim(data)[1],ncol=nocc-1)
   age.mat=t(t(age.mat)+cumsum(c(0,time.intervals[1:(nocc-2)])))
   age.mat=data.frame(age.mat)
   names(age.mat)=paste("NestAge",1:(nocc-1),sep="")
   return(age.mat)
}


#
#  *******************  END OF INTERNAL FUNCTIONS    *********************************
# Test to make sure that all rows of design data are there (no more deletion) and make sure they
# are ordered
  for(i in 1:(length(ddl)-1))
  {
	  if(max(ddl[[i]]$par.index) != nrow(ddl[[i]])) 
	  {
		  warning(paste("\nMissing rows in design dataframe for parameter",names(ddl)[i],
				  "\n Deleting rows from design data is still allowed but see warning in help for make.design.data\n"))
	  }
	 if(any(ddl[[i]]$par.index != sort(ddl[[i]]$par.index))) 
	 {
		 stop(paste("\nRows of design dataframe for parameter",names(ddl)[i],
				 "are out of order.\nThey should be in order of par.index.\n"))
 	 }
  }
#
# Outfile is assigned a temporary name
#
  tempfilename=tempfile("markxxx",tmpdir=getwd(),fileext=".tmp")
  outfile=file(tempfilename,open="wt")
#
# Check validity of parameter types, if any given
#
  if(!valid.parameters(data$model,parameters))stop()
#
# Next check validity of fields defined in each parameter list.
#
  if(length(parameters)>0)
  for(i in 1:length(parameters))
     for(j in 1:length(names(parameters[[i]])))
        if(!(names(parameters[[i]])[j]%in%c("fixed","formula","link","share","remove.intercept","default")))
        {
           message("\nInvalid model specification for parameter ",names(parameters)[i],".\nUnrecognized element ",names(parameters[[i]])[j])
           stop()
        }
#
# Initialize some variables
#
  ch=data$data$ch
  mixtures=data$mixtures
  nocc=data$nocc
  nocc.secondary=data$nocc.secondary
  nstrata=data$nstrata
  number.of.groups=dim(data$freq)[2]
  par.list=setup.parameters(data$model,check=TRUE)
  parameters=setup.parameters(data$model,parameters,nocc,number.of.groups=number.of.groups)
  parameters=parameters[par.list]
  temp.rev=data$reverse
  data$reverse=FALSE
  full.ddl=make.design.data(data,parameters=ddl$pimtypes)
  data$reverse=temp.rev
  parameters=parameters[names(parameters)%in%names(full.ddl)]
  for(j in names(parameters))
  {
     parameters[[j]]$pim.type=ddl$pimtypes[[j]]$pim.type
     if(!is.null(ddl$pimtypes[[j]]$subtract.stratum))
        parameters[[j]]$subtract.stratum=ddl$pimtypes[[j]]$subtract.stratum
  }
  for(i in 1:length(parameters))
  {
#
#     For parameters that can be possibly shared, see if they are not shared and if not then create
#     default formula if one not specified; also use link from dominant parameter
#
	  if(!is.null(parameters[[i]]$share))
	  {
		  if(!parameters[[i]]$share)
          {
		      shared_par=parameters[[i]]$pair
	          if(is.null(parameters[[shared_par]]$formula))parameters[[shared_par]]$formula=~1
	      }else
		  {
			  shared_par=parameters[[i]]$pair
			  parameters[[shared_par]]$link=parameters[[i]]$link	  
		  }
	  }
#
#     Test validity of link functions
# 
     if(!(parameters[[i]]$link %in% c("identity","log","logit","mlogit","loglog","cloglog","sin")))
     {
        message("\nInvalid link value ",parameters[[i]]$link)
        stop()
     }
  }
  param.names=sub("DoublePrime","''",names(parameters))
  param.names=sub("Prime","'",param.names)
  model.list= setup.model(data$model,0)
  etype=model.list$etype
#	
# Output data portion of MARK input file:
#   proc title
#   proc chmatrix
#
  if(etype=="Nest")
  {
     zz=subset(data$data,select=c("FirstFound","LastPresent","LastChecked","Fate"))
     zz=cbind(zz,rowSums(data$freq))
     if(!is.null(data$data$AgeDay1))
        data$data=cbind(data$data,create.agenest.var(data$data,"AgeDay1",data$time.intervals))
  }
  else
  {
     zz=as.data.frame(ch)
     zz=cbind(zz,data$freq)
  }
#
# p&c in closed models, gammas in robust models, and p1-p2 in MS Occupancy are handled differently
# to allow shared parameters.  If there is no c$formula, gammaDoublePrime$formula, or p2$formula then the design data
# are appended to the shared parameters (p for c and gammaPrime for gammaDoublePrime)
# In the design data for p a covariate "c" is added to the recapture parameters and for
# a covariate emigrate for the gammaDoublePrime parameters.
# 
    for(i in 1:length(parameters))
   {
#
#     For parameters that can be possibly shared, if they are shared, pool design data as long as dimensions match
#
	   if(!is.null(parameters[[i]]$share)&&parameters[[i]]$share)
       {
	       shared_par=parameters[[i]]$pair
	       dim1=dim(ddl[[names(parameters)[i]]])
		   dim2=dim(ddl[[shared_par]])
           if(dim1[2]==dim2[2])
		   {
			   rn1=as.numeric(row.names(ddl[[names(parameters)[i]]]))
			   rn2=as.numeric(row.names(ddl[[shared_par]]))+nrow(full.ddl[[names(parameters)[i]]])
			   ddl[[names(parameters)[i]]]=rbind(ddl[[names(parameters)[i]]],ddl[[shared_par]])
               ddl[[names(parameters)[i]]][shared_par]=c(rep(0,dim1[1]),rep(1,dim2[1]))
		       row.names(ddl[[names(parameters)[i]]])=c(rn1,rn2)
		   } else
		   {
			   message(paste("Error: for a shared ",paste(names(parameters)[i],shared_par,sep="&"),
				" model, their design data columns must match\n. If you add design data to one it must also be added to the other.\n"))
			   message(paste("Columns of",names(parameters)[i]," : ",names(ddl[i]),"\n"))
			   message(paste("Columns of",shared_par,": ",names(ddl[[shared_par]]),"\n"))
			   stop("Function terminated\n") 
		   }
	   }
   }
#
# For each parameter type determine which values in the formula are covariates that need to be
# added to design data and put in data portion of input file.
#
  xcov=list()
  covariates=NULL
  time.dependent=list()
  session.dependent=list()
  for(i in 1:length(parameters))
  {
     if(!is.null(parameters[[i]]$formula))
     {
#
#    First get the variables in the formula.  Identify those not in the design data and make sure
#    that they are in the data. If there are any, add the covariate to the covariate list and add a column 
#    to the design data.
#
#    Note parx is the name of the ith parameter.  ddl is always constructed in the same order for each model
#    but the order of the parameters in the model specification may be different, therefore the indexing into
#    ddl must be done by name rather than position from the parameters specification.
#
      parx=names(parameters)[i]
      vars=all.vars(parameters[[i]]$formula)
      termslist=colnames(attr(terms(parameters[[i]]$formula),"factors"))
      xcov[[parx]]=vars[!(vars%in%names(ddl[[parx]]))]
      time.dependent[[parx]]=rep(FALSE,length(xcov[[parx]]))
      session.dependent[[parx]]=rep(FALSE,length(xcov[[parx]]))
      if(any(!(vars%in%names(ddl[[parx]]))))
         for(j in 1:length(xcov[[parx]]))
         {
            if(!(xcov[[parx]][j]%in%names(data$data)))
            {
              if(!is.null(full.ddl[[parx]]$session))
              {
                 cov.bytime=unique(paste(xcov[[parx]][j],as.character(ddl[[parx]]$session),as.character(ddl[[parx]]$time),sep=""))
                 if(any(!cov.bytime%in%names(data$data)))
                 {
                    session.dependent[[parx]][j]=TRUE                 
                    cov.bytime=unique(paste(xcov[[parx]][j],as.character(ddl[[parx]]$session),sep=""))
                 }
              }                 
              else
                 cov.bytime=unique(paste(xcov[[parx]][j],as.character(ddl[[parx]]$time),sep=""))
              if(all(!cov.bytime%in%names(data$data)))
                 stop("\nError: Variable ",xcov[[parx]][j]," used in formula is not defined in data\n")
              else
              {
                 time.dependent[[parx]][j]=TRUE
                 if(any(!cov.bytime%in%names(data$data)))
                 {
                     message(paste("\nThe following covariates are missing:",cov.bytime[!cov.bytime%in%names(data$data)],collapse=""))
                     message("\nIf any are used in the resulting model it will fail\n")
                     cov.bytime=cov.bytime[cov.bytime%in%names(data$data)]
                 }
              }
            }
            savenames=names(ddl[[parx]])
            ddl[[parx]]=cbind(ddl[[parx]],rep(1,dim(ddl[[parx]])[1]))
            names(ddl[[parx]])=c(savenames,xcov[[parx]][j])
            if(!time.dependent[[parx]][j])
               covariates=c(covariates,xcov[[parx]][j])
            else
               covariates=c(covariates,cov.bytime)
         }
     }
  }
# 
#  Unless this is nest data, aggregate data
#
# Fix 9 Nov 2013; create unique covariate names before selecting data
  covariates=unique(covariates)
  if(!is.null(covariates))
	    zzd=data.frame(cbind(zz,data$data[,covariates]))
  else
	    zzd=zz
  if(etype!="Nest" & accumulate)
  {
	  pasted.data=apply(zzd, 1, paste, collapse = "")
	  ng=ncol(data$freq)
	  if(ng>1)
		  freq=t(sapply(split(data$freq,pasted.data ),colSums))
	  else
		  freq=sapply(split(data$freq, pasted.data),sum)
	  zzd=zzd[order(pasted.data),]
	  zzd=zzd[!duplicated(pasted.data[order(pasted.data)]),]
	  if(ng>1)
	  {
		  if(nrow(freq)!=nrow(zzd))
			  stop("problem with accumulating data. Set accumulate=FALSE and contact maintainer")
	  }else
	  if(length(freq)!=nrow(zzd))
		  stop("problem with accumulating data. Set accumulate=FALSE and contact maintainer")
	  zzd[,2:(ng+1)]=freq
	  zz=zzd[,1:(ng+1)]
  }
#
# Output title and list of covariates
# 11 Jan 06; Added code for multistratum - nstrata and strata labels
#

  if(is.null(nocc.secondary))
     string=paste("proc title ",title,";\nproc chmatrix occasions=",nocc," groups=",number.of.groups," etype=",etype)
  else
     string=paste("proc title ",title,";\nproc chmatrix occasions=",sum(nocc.secondary)," groups=",number.of.groups," etype=",etype)
  if(model.list$strata)string=paste(string," strata=",data$nstrata,sep="")
  if(!is.null(covariates))
  {
	 covar10=covariates[duplicated((substr(covariates,1,10)))]
	 if(length(covar10)>0) stop(paste("\nFollowing covariates are duplicates of another covariate within the first 10 characters\n",paste(covar10,collapse=", ")))
     string=paste(string," icovar = ",length(covariates))
	 if(!is.null(icvalues))
	 {
		 if(length(covariates)!=length(icvalues))
			 stop("\nMismatch between length of individual covariates and covariate values (icvalues)\n")
	     if(!is.numeric(icvalues)) 
			 stop("\nicvalues must be numeric\n")
	 }
  }
  if(mixtures!=1)
      string=paste(string," mixtures =",mixtures)
  if(data$model=="MultScalOcc")
	  string=paste(string," mixtures =",length(levels(ddl$p$time)))
  time.int=data$time.intervals
  if(!is.null(data$reverse) &&(data$reverse | data$model=="MultScalOcc")) time.int[time.int==0]=1
  string=paste(string," ICMeans NoHist hist=",nrow(zz),
           ";\n time interval ",paste(time.int,collapse=" "),";\n")
  if(model.list$strata)string=paste(string,"strata=",paste(data$strata.labels[1:data$nstrata],collapse=" "),";\n",sep="")
  if(!is.null(covariates))
  {
     string=paste(string,"icovariates ",paste(covariates,collapse=" "),";")
     any.factors=sapply(data$data[,covariates,drop=FALSE],is.factor)
     if(any(any.factors))
        stop(paste("The following individual covariates are not allowed because they are factor variables: ",paste(names(data$data[,covariates,drop=FALSE])[any.factors],collapse=",")))
     else
     {
        any.na=apply(data$data[,covariates,drop=FALSE],2,function(x) any(is.na(x)))
        if(any(any.na))
           stop(paste("The following individual covariates are not allowed because they contain NA: ",paste(names(data$data[,covariates,drop=FALSE])[any.na],collapse=",")))
        else
            zz=zzd
     }
  }
  write(strwrap(string,100,prefix=" "),file=outfile,append=TRUE)
#
#  Output group labels
#
  if(is.null(names(data$freq)))
     group.labels=paste("Group",1:number.of.groups)
  else
     group.labels=names(data$freq) 
  for(j in 1:number.of.groups)
  {
      string=paste("glabel(",j,")=",group.labels[j],";",sep="")
      write(string,file=outfile,append=TRUE)
  }
#
# This outputs capture history, frequency and any covariates
#
  if(etype=="Nest")
  {
     for(i in 1:number.of.groups)
     {
        string=paste("Nest Survival Group =", i, ";")
        write(string,file=outfile,append=TRUE)
        if(number.of.groups>1)
            write.table(zz[data$data$group==i,],file=outfile,eol=";\n",sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
        else
            write.table(zz,file=outfile,eol=";\n",sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
     }
  }
  else{
		if(!wrap)
		     write.table(zz,file=outfile,eol=";\n",sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
        else 
		   apply(zz,1,function(x) write(strwrap(paste(paste(x,collapse=" "),";",sep=""),80,prefix=" "),file=outfile,append=TRUE))
   }
#
# Output counts section of Mark-resight models if appropriate
#
  if(!is.null(data$counts))
	  for(i in 1:number.of.groups)
	  {
		  for(j in names(data$counts))
		  {
			  write("",file=outfile,append=TRUE)
			  string=paste(j,"Group =", i, ";")
		      write(string,file=outfile,append=TRUE)
		      if(number.of.groups==1)
		          string=paste(paste(data$counts[[j]],collapse=" "),";",sep="")
		      else
			      string=paste(paste(data$counts[[j]][i,],collapse=" "),";",sep="")
 		      write(string,file=outfile,append=TRUE)
		  }
	  }  
#
# First create model name using each parameter unless a model name was given as an argument;
#
  if(is.null(model.name))
  {
    model.name=""
    for(i in 1:length(parameters))
    {
       model.name=paste(model.name,param.names[i],"(",paste(parameters[[i]]$formula,collapse=""),sep="")
       model.name=paste(model.name,")",sep="")
     }
  }
#
# Next determine if a single link was used or differing links
#
  link=parameters[[1]]$link
  for(i in 1:length(parameters))
     if(link!=parameters[[i]]$link)link="Parm-Specific"
  if(!is.null(input.links) | parm.specific)link="Parm-Specific"
#
# Next output proc Estimate 
#
  if(!hessian)
     string=paste(paste("proc estimate link=",spell(link),sep=""),"NOLOOP varest=2ndPart ",options," ;\nmodel={",substr(model.name,1,160),"};")
  else
     string=paste(paste("proc estimate link=",spell(link),sep=""),"NOLOOP varest=Hessian ",options," ;\nmodel={",substr(model.name,1,160),"};")
 write(string,file=outfile,append=TRUE)
#
# Next compute PIMS for each parameter in the model - these are the all different PIMS
# 11 Jan 06 code added for multistrata models
#
  pim=list()
  npar=1
  for(i in 1:length(parameters))
  {
     pim[[i]]=list()
     k=0
     for(j in 1:number.of.groups)
     {
	   if(is.null(parameters[[i]]$bystratum)||!parameters[[i]]$bystratum)
         xstrata=1
	   else
		 xstrata=unique(ddl[[i]]$stratum)
	   for (jj in xstrata)
	   {
          other.strata=1
          if(!is.null(parameters[[i]]$tostrata))
			  other.strata=unique(ddl[[i]]$tostratum[ddl[[i]]$stratum==jj])
          for(to.stratum in other.strata)
          {
               if(model.list$robust && parameters[[i]]$secondary)
			   {
				   multi.session=TRUE
                   num.sessions=nocc
			   }
               else
			   {
				   num.sessions=1
				   multi.session=FALSE
			   }
               for (l in 1:num.sessions)
               {
                  k=k+1
                  pim[[i]][[k]]=list()
				  if(data$model=="RDMSOccRepro" & names(parameters)[i]=="Phi0")
				  {
					  pim[[i]][[k]]$pim=matrix(1:length(data$strata.labels),nrow=1)   
				  } else	 
				  if(!multi.session)
					 pim[[i]][[k]]$pim=create.pim(nocc,parameters[[i]],npar,mixtures)
				  else
                  {
                     if(is.na(parameters[[i]]$num))
                     {
                         parameterx=parameters[[i]]
                         parameterx$num=0
                         pim[[i]][[k]]$pim=create.pim(1,parameterx,npar,mixtures)
                     }
                     else
                         pim[[i]][[k]]$pim=create.pim(nocc.secondary[l],parameters[[i]],npar,mixtures)
                     pim[[i]][[k]]$session=l
					 pim[[i]][[k]]$session.label=levels(ddl[[i]]$session)[l]
                  }
                  pim[[i]][[k]]$group=j
                  if(length(data$strata.labels)>0 && !is.null(parameters[[i]]$bystratum) && parameters[[i]]$bystratum) pim[[i]][[k]]$stratum=jj
                  if(!is.null(parameters[[i]]$tostrata)) pim[[i]][[k]]$tostratum=to.stratum
                  npar=max(pim[[i]][[k]]$pim)+1
               }
           }
       }
     }
  }
  npar=npar-1
  names(pim)=names(parameters)
#
# If there are fixed parameters output the text here
#
  num.fixed=0
  fixed=NULL
  fixedvalues=NULL
  for(i in 1:length(parameters))
  {
     parx=names(parameters)[i]
#
#    Add any default fixed values
#
     fixlist=NULL
     fixvalues=NULL
     fix.indices=NULL
     if(default.fixed)
     {
        rn=row.names(full.ddl[[parx]][!row.names(full.ddl[[parx]])%in%row.names(ddl[[parx]]),])
        if(length(rn)>0)
        {
           fixvalues=rep(parameters[[i]]$default,length(rn))
           fixlist=as.numeric(rn)
        }
     }
#
#    Add any values specified with fix column in ddl
#
	 if(!is.null(ddl[[parx]]$fix))
	 {
		 fixvalues=c(fixvalues,ddl[[parx]]$fix[!is.na(ddl[[parx]]$fix)])
		 fixlist=c(fixlist,as.numeric(row.names(ddl[[parx]][!is.na(ddl[[parx]]$fix),])))
	 }
     if(!is.null(parameters[[i]]$fixed)|!is.null(fixlist))
     {
#
#      All values of this parameter type are fixed at one value
#
       if(length(parameters[[i]]$fixed)==1) 
       {
          for(j in 1:length(pim[[i]]))
          {
             fixlist=unique(as.vector(pim[[i]][[j]]$pim))
             fixlist=fixlist[fixlist>0]
             if(is.null(fixedvalues))
                 fixedvalues=data.frame(index=fixlist,value=rep(parameters[[i]]$fixed,length(fixlist)))
             else
                 fixedvalues=rbind(fixedvalues,data.frame(index=fixlist,value=rep(parameters[[i]]$fixed,length(fixlist))))
             for(k in 1:length(fixlist))
             {
                num.fixed=num.fixed+1
                fixed=c(fixed,paste("parm(",fixlist[k],")=",parameters[[i]]$fixed,sep=""))
             }
          }
        }else
#
#       The parameters with indices in the first list element are given specified value(s)
#
        if(is.list(parameters[[i]]$fixed)|!is.null(fixlist))
        {
             if("index"%in%names(parameters[[i]]$fixed))
                fix.indices=parameters[[i]]$fixed$index
             else
                if(is.list(parameters[[i]]$fixed)&!any(names(parameters[[i]]$fixed)%in%c("time","age","cohort","group")))
                   stop(paste("\nUnrecognized structure for fixed parameters =",parameters[[i]]$fixed))
                else
                    if(!is.null(parameters[[i]]$fixed[["time"]]))
                    {
                        times=parameters[[i]]$fixed[["time"]]
                        fix.indices=as.numeric(row.names(full.ddl[[parx]][full.ddl[[parx]]$time%in%times,]))
                    }
                    else
                    if(!is.null(parameters[[i]]$fixed[["age"]]))
                    {
                        ages=parameters[[i]]$fixed[["age"]]
                        fix.indices=as.numeric(row.names(full.ddl[[parx]][full.ddl[[parx]]$age%in%ages,]))
                    }
		                else
                    if(!is.null(parameters[[i]]$fixed[["cohort"]]))
                    {
                        cohorts=parameters[[i]]$fixed[["cohort"]]
                        fix.indices=as.numeric(row.names(full.ddl[[parx]][full.ddl[[parx]]$cohort%in%cohorts,]))
                    }
                    else
                    if(!is.null(parameters[[i]]$fixed[["group"]]))
                    {
                        groups=parameters[[i]]$fixed[["group"]]
                        fix.indices=as.numeric(row.names(full.ddl[[parx]][full.ddl[[parx]]$group%in%groups,]))
                    }
             if(!is.null(fix.indices))fixlist=c(fixlist,fix.indices)
             if(length(parameters[[i]]$fixed$value)==1)
                 fixvalues=c(fixvalues,rep(parameters[[i]]$fixed$value,length(fix.indices)))
             else
             {
                 fixvalues=c(fixvalues,parameters[[i]]$fixed$value)
                 if(length(fixlist)!=length(fixvalues))
                    stop(paste("\nLengths of indices and values do not match for fixed parameters for",names(parameters)[i],"\n"))
             }
			 # check for duplicates and use latter values
             if(any(duplicated(fixlist)))
			 {
				 message(paste("\nSome indices for fixed parameters were duplicated for ",parx,"; using latter values\n"))
				 uniqIndices=which(!duplicated(rev(fixlist)))
				 fixlist=rev(fixlist)[uniqIndices]
				 fixvalues=rev(fixvalues)[uniqIndices]
			 }
             # assign all.different indices by adding first pim index-1
             fixlist=fixlist+ pim[[i]][[1]]$pim[1,1]-1
             for(k in 1:length(fixlist))
             {
                num.fixed=num.fixed+1
                fixed=c(fixed,paste("parm(",fixlist[k],")=",fixvalues[k],sep=""))
             }
             if(is.null(fixedvalues))
                fixedvalues=data.frame(index=fixlist,value=fixvalues)
             else
                fixedvalues=rbind(fixedvalues,data.frame(index=fixlist,value=fixvalues))
        }
     }
  }
  if(num.fixed>0)
  {
     string=paste("fixed =",num.fixed,";",sep="")
     write(string,file=outfile,append=TRUE)
     write(paste(fixed,";"),file=outfile,append=TRUE)
  }
#
# Unless model will be simplified, output PIMS for each parameter in the model
#
  if(!simplify)
  {
     for(i in 1:length(parameters))
     {
        for(j in 1:length(pim[[i]]))
        {
           ncol=dim(pim[[i]][[j]]$pim)[2]
           string=pim.header(pim[[i]][[j]]$group,param.names[i],parameters[[i]],
                   ncol,pim[[i]][[j]]$stratum,pim[[i]][[j]]$tostratum,
                   data$strata.labels,mixtures,pim[[i]][[j]]$session,parameters[[i]]$socc)
           write(string,outfile,append=TRUE)
           print.pim( pim[[i]][[j]]$pim,outfile)
        }
     }
  }
#
# Create design matrix for each parameter type
#
  design.matrix=list()
  for(i in 1:length(parameters))
  {
  if(is.null(parameters[[i]]$formula))
  {
     design.matrix[[i]]=list()
  } else
  {    
#  Next, if share, combine full ddl
	fullddl=full.ddl[[names(parameters)[i]]]
	if(!is.null(parameters[[i]]$share)&&parameters[[i]]$share)
	{	  
		  shared_par=parameters[[i]]$pair
		  dim1=dim(fullddl)
		  dim2=dim(full.ddl[[shared_par]])
		  fullddl=rbind(fullddl,full.ddl[[shared_par]])
		  fullddl[shared_par]=c(rep(0,dim1[1]),rep(1,dim2[1]))
		  row.names(fullddl)=1:dim(fullddl)[1]
	  } 
#
#    Calculate number of parameters for this type
#
     parx=names(parameters)[i]
     npar=dim(ddl[[parx]])[1]
#
#       Compute design matrix with model.matrix
#         31 Jan 06; made change to allow for NA or missing design data
#
        dm=model.matrix(parameters[[parx]]$formula,ddl[[parx]])
#
#       In cases with nested interactions it is necessary to remove the intercept
#       to avoid over-parameterizing the model; this is user-specified
#
        if(!is.null(parameters[[parx]]$remove.intercept)&&parameters[[parx]]$remove.intercept)
        {
            intercept.column=(1:dim(dm)[2])[colSums(dm)==dim(dm)[1]]
            if(length(intercept.column)==0)
               stop("\nIntercept column not found.  Do not use ~-1 with remove.intercept\n")
            else
            {
               if(length(intercept.column)==1)
                  dm=dm[,-intercept.column]
            }
        }
        maxpar=dim(fullddl)[1]
#		Create a complete design matrix using full ddl
        design.matrix[[i]]=matrix(0,ncol=dim(dm)[2],nrow=maxpar)
#       Using row numbers fill in the dm for rows design data; this handles deleted design data
		if(dim(design.matrix[[i]][as.numeric(row.names(ddl[[parx]])),,drop=FALSE])[1]==dim(dm)[1])
           design.matrix[[i]][as.numeric(row.names(ddl[[parx]])),]=dm
        else
           stop(paste("\nProblem with design data. It appears that there are NA values in one or more variables in design data for ",parx,"\nMake sure any binned factor completely spans range of data\n",sep=""))
        colnames(design.matrix[[i]])=colnames(dm)
#
#       It appears that model.matrix can add unneeded columns to the design matrices
#       It can add interactions that are not relevant.  The results are columns in the design
#       matrix that are all zero.  These are stripped out here.
#
        col.sums=apply(design.matrix[[i]],2,function(x) sum(abs(x)))
        design.matrix[[i]]=design.matrix[[i]][,col.sums!=0,drop=FALSE]
#
#       Next substitute variable names for covariates into design matrix 
#
        if(length(xcov[[parx]])!=0)
          for(j in 1:length(xcov[[parx]]))
          {
             which.cols=NULL
             cnames=colnames(design.matrix[[i]])
#
#            fix for v1.6.2 needed to delete [ or ] or ( or ) or , from the names
#            to use the all.vars command.  These are created in using the cut command
#            on a numeric variable.
#
             cnames=sub("\\(","",cnames)
             cnames=sub("\\)","",cnames)
             cnames=sub("\\[","",cnames)
             cnames=sub("\\]","",cnames)
             cnames=sub(",","",cnames)
             for(k in 1:dim(design.matrix[[i]])[2])
                if(xcov[[parx]][j]%in%all.vars(formula(paste("~",cnames[k],sep=""))))
                     which.cols=c(which.cols,k)
             if(length(which.cols)>0)
             for(k in which.cols)
               if(all(design.matrix[[i]][,k]==1 | design.matrix[[i]][,k]==0))
                  if(time.dependent[[parx]][j])
                     if(!is.null(fullddl$session))
                        if(session.dependent[[parx]][j])
                            design.matrix[[i]][,k][design.matrix[[i]][,k]==1]=
                                paste(xcov[[parx]][j],as.character(fullddl$session[design.matrix[[i]][,k]==1]),sep="")
                        else
                           design.matrix[[i]][,k][design.matrix[[i]][,k]==1]=
                                paste(xcov[[parx]][j],as.character(fullddl$session[design.matrix[[i]][,k]==1]),
                                  as.character(fullddl$time[design.matrix[[i]][,k]==1]),sep="")
                     else   
                        design.matrix[[i]][,k][design.matrix[[i]][,k]==1]=paste(xcov[[parx]][j],as.character(fullddl$time[design.matrix[[i]][,k]==1]),sep="")
                  else
                     design.matrix[[i]][,k][design.matrix[[i]][,k]==1]=xcov[[parx]][j]
               else
                  if(time.dependent[[parx]][j])
                  {
                     if(!is.null(fullddl$session))
                     {
                        if(session.dependent[[parx]][j])
                           design.matrix[[i]][,k][design.matrix[[i]][,k]!=0]=
                              paste("product(",paste(xcov[[parx]][j],as.character(fullddl$session[design.matrix[[i]][,k]!=0]),sep=""),
                                       ",",design.matrix[[i]][,k][design.matrix[[i]][,k]!=0],")",sep="")
                        else
                           design.matrix[[i]][,k][design.matrix[[i]][,k]!=0]=
                              paste("product(",paste(xcov[[parx]][j],as.character(fullddl$session[design.matrix[[i]][,k]!=0]),
                                      as.character(fullddl$time[design.matrix[[i]][,k]!=0]),sep=""),
                                       ",",design.matrix[[i]][,k][design.matrix[[i]][,k]!=0],")",sep="")                     
                     }
                     else
                        design.matrix[[i]][,k][design.matrix[[i]][,k]!=0]=
                           paste("product(",paste(xcov[[parx]][j],as.character(fullddl$time[design.matrix[[i]][,k]!=0]),sep=""),
                                    ",",design.matrix[[i]][,k][design.matrix[[i]][,k]!=0],")",sep="")
                  }
                  else
                     design.matrix[[i]][,k][design.matrix[[i]][,k]!=0]=
                        paste("product(",xcov[[parx]][j],",",design.matrix[[i]][,k][design.matrix[[i]][,k]!=0],")",sep="")
          }

        row.names(design.matrix[[i]])=NULL
     design.matrix[[i]]=as.data.frame(design.matrix[[i]])
     if(parameters[[i]]$formula=="~1")
        names(design.matrix[[i]])[1]="(Intercept)"
     names(design.matrix[[i]])=paste(names(parameters)[i],names(design.matrix[[i]]),sep=":")
  } 
  }
  names(design.matrix)=names(parameters)
#
# Merge to create a single design matrix
#
  complete.design.matrix=NULL
  nrows=0
  lastpim=length( pim[[length(parameters)]])
  lastindex=sum(sapply(full.ddl[1:length(parameters)],nrow))
#  lastindex=max(pim[[length(parameters)]][[lastpim]]$pim)
  for(i in 1:length(parameters))
  {
	 # parameters with NULL formula have been merged with a shared parameter
     if(!is.null(parameters[[i]]$formula))
     {    
        mat=NULL
		pair=parameters[[i]]$pair
		if(!is.null(pair) && pair !="" && !is.null(parameters[[pair]]$share) && parameters[[pair]]$share)
		{
			minrow=pim[[names(parameters)[i]]][[1]]$pim[1,1]
			maxrow=max(pim[[names(parameters)[i]]][[length(pim[[names(parameters)[i]]])]]$pim)
			if(minrow>1)
				mat=matrix("0",ncol=dim(design.matrix[[i]])[2],nrow=minrow-1)
			mat=rbind(mat,as.matrix(design.matrix[[i]]))
			if(i<length(parameters))
				mat=rbind(mat,matrix("0",ncol=dim(design.matrix[[i]])[2],nrow=lastindex-maxrow ))
		} else
        {    
           if(i>1)
              mat=matrix("0",ncol=dim(design.matrix[[i]])[2],nrow=nrows)
           mat=rbind(mat,as.matrix(design.matrix[[i]]))
           nrows=dim(mat)[1]
           if(i<length(parameters))
              mat=rbind(mat,matrix("0",ncol=dim(design.matrix[[i]])[2],nrow=lastindex-nrows ))
        }
        names(mat)=names(design.matrix[[i]])
        complete.design.matrix=cbind(complete.design.matrix,mat)
     }
  }
  row.names(complete.design.matrix)=1:dim(complete.design.matrix)[1]
  complete.design.matrix=as.data.frame(complete.design.matrix)

#
#  If there are any initial values, output them to the MARK input file
#  after making sure that the vector length matches the number of parameters  
#
   if(!is.null(initial))
   {
#
#     If a vector of values was given check to make sure it is of the correct length and then output
#
      if(is.vector(initial))
      {
         if(length(initial)==dim(complete.design.matrix)[2])
            initial.values=initial
         else
         {
            if(length(initial)==1)
               initial.values=rep(initial,dim(complete.design.matrix)[2])
            else
            {
                if(length(names(initial))==0)
				{
				  message("\nLength of initial vector doesn't match design matrix: ",ncol(complete.design.matrix)," \n")
				  print(colnames(complete.design.matrix))
                  stop()
			    }else
				{
					beta.index=match(names(complete.design.matrix),names(initial))
					initial.values=rep(0,dim(complete.design.matrix)[2])
					initial.values[!is.na(beta.index)]=initial[beta.index[!is.na(beta.index)]]				
				}
            }
         }
      } 
      else
#
#     If it was a MARK object; check to make sure it has output and then use the betas from the other object
#     2 May 06 jll; use names instead of values in design matrix to match initial values
#
      {
         if(class(initial)[1]=="mark")
         {
            initial=load.model(initial)
            if(!is.null(initial$output))
            {
               beta.index=match(names(complete.design.matrix),colnames(initial$design.matrix))
               initial.values=rep(0,dim(complete.design.matrix)[2])
               initial.values[!is.na(beta.index)]=initial$results$beta$estimate[beta.index[!is.na(beta.index)]]
            }
         }
      }
      if(simplify)
		  string=paste("XXXinitialXXX ",paste(initial.values,collapse=" "),";")
	  else
		  string=paste("initial ",paste(initial.values,collapse=" "),";")
	  write(string,file=outfile,append=TRUE)
   }
#
#  If model will not be simplified, output design matrix to the MARK input file
#
  if(!simplify)
  {
    string=paste("design matrix constraints=",dim(complete.design.matrix)[1], " covariates=",dim(complete.design.matrix)[2],";",sep="")
    write(string,file=outfile,append=TRUE)
    write.table(complete.design.matrix,file=outfile,eol=";\n",sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
  }
#
#  If parm-specific links, output them here
# 
  mlogit.list=list(structure=NULL,ncol=1)
  if(is.null(input.links))
  {
     max.logit.number=0
     if(link=="Parm-Specific")
     {
        string=NULL
        for(i in 1:length(parameters))
        {
           parx=names(parameters)[i]
           if(parameters[[i]]$link=="mlogit"|parameters[[i]]$link=="MLogit")
           {
             if(parx=="Psi")
             {
                 logit.numbers = max.logit.number+1:(nrow(full.ddl[[parx]])/(nstrata*(nstrata-1)*number.of.groups))
                 logits.per.group=nstrata*length(logit.numbers)
                 for (k in 1:number.of.groups)
                 {
                    if(k>1)logit.numbers=logit.numbers+logits.per.group
                    for (j in 1:nstrata)	 
                      string=c(string,paste("mlogit(",rep(logit.numbers+(j-1)*length(logit.numbers),(nstrata-1)),")",sep="")) 
			     }
			     max.logit.number=max.logit.number+logits.per.group*number.of.groups
              }
              else
              {
                 if(parx%in% c("pent","alpha"))
                 {
                     nsets=length(pim[[parx]])
                     for (kk in 1:nsets)
                     {
                        x.indices=as.vector(t(pim[[parx]][[kk]]$pim))
                        x.indices=x.indices[x.indices!=0]
                        max.logit.number=max.logit.number+1
                       string=c(string,paste("mlogit(",rep(max.logit.number,length(x.indices)),")",sep=""))
                     }
                 }else
                 {
				    if(parx %in% c("pi","Omega"))
				    {
					    for (kk in 1:number.of.groups)
					    {
					   	     logit.numbers=max.logit.number+rep(1:(nrow(full.ddl[[parx]])/(number.of.groups*(nstrata-1))),nstrata-1)
						     max.logit.number=max(logit.numbers)
						     string=c(string,paste("mlogit(",logit.numbers,")",sep=""))
					    }			 				 
				     } else
                         stop(paste("Mlogit link not allowed with parameter",parx))
			     }
               }
            } else
            {
              xstring=rep(spell(parameters[[i]]$link),dim(full.ddl[[parx]])[1])
              string=c(string,xstring)
            }
        }
   	    write(paste("links=",length(string),";",sep=""),file=outfile,append=TRUE)
        links=string
        string=paste(string,";")
        write(string,file=outfile,append=TRUE)
     } else
     {
		 links=link
     }
  }	else
  {
	  string=paste(input.links,collapse=",")
	  write(paste("links=",length(string),";",sep=""),file=outfile,append=TRUE)
	  links=input.links
	  string=paste(string,";")
	  write(string,file=outfile,append=TRUE)
  }  
#
# write out labels for design matrix columns
#
  for(i in 1:dim(complete.design.matrix)[2])
  {
     string=paste("blabel(",i,")=",colnames(complete.design.matrix)[i],";",sep="")
     write(string,file=outfile,append=TRUE)
  }
#
# write out labels for real parameters
#
  labstring=NULL
  rnames=NULL
  ipos=0
  for(i in 1:length(parameters))
  {
      parx=names(parameters)[i]
      plimit=dim(full.ddl[[parx]])[1]
      stratum.strings=rep("",plimit)
      if(!is.null(full.ddl[[parx]]$stratum)) stratum.strings=paste(" s",full.ddl[[parx]]$stratum,sep="")
      if(!is.null(full.ddl[[parx]]$tostratum)) stratum.strings=paste(stratum.strings," to",full.ddl[[parx]]$tostratum,sep="")
      strings=paste(param.names[i],stratum.strings," g",full.ddl[[parx]]$group,sep="")
      if(data$reverse)
	  {
		  if(!is.null(full.ddl[[parx]]$cohort))strings=paste(strings," c",full.ddl[[parx]]$cohort,sep="")
		  if(!is.null(full.ddl[[parx]]$occ.cohort))strings=paste(strings," c",full.ddl[[parx]]$occ.cohort,sep="")
	  } else {
		  if(!is.null(full.ddl[[parx]]$cohort))
			  strings=paste(strings," c",full.ddl[[parx]]$cohort,sep="")
		  else
		  if(!is.null(full.ddl[[parx]]$occ.cohort))strings=paste(strings," c",full.ddl[[parx]]$occ.cohort,sep="")
	  }  
	  if(!is.null(full.ddl[[parx]]$age))strings=paste(strings," a",full.ddl[[parx]]$age,sep="")
	  if("occ"%in%names(full.ddl[[parx]]))strings=paste(strings," o",full.ddl[[parx]]$occ,sep="")
	  if(model.list$robust && parameters[[parx]]$secondary)
         strings=paste(strings," s",full.ddl[[parx]]$session,sep="")
      if(!is.null(full.ddl[[parx]]$time))strings=paste(strings," t",full.ddl[[parx]]$time,sep="")
      if(mixtures >1 && !is.null(parameters[[i]]$mix) &&parameters[[i]]$mix)
         strings=paste(strings," m",full.ddl[[parx]]$mixture,sep="")
      rnames=c(rnames,strings)
      if(!simplify)
      {
         strings=paste("rlabel(",ipos+1:plimit,")=",strings,";",sep="")
         labstring=c(labstring,strings)
         ipos=ipos+plimit
      }
  }
  if(any(duplicated(rnames))) stop("Contact package maintainer. Following row names are duplicated:",paste(rnames[duplicated(rnames)],collapse=" "))
  if(!simplify) write(labstring,file=outfile,append=TRUE)
  row.names(complete.design.matrix)=rnames
#
# Write out Proc stop statement
#
  write("proc stop;",file=outfile,append=TRUE)
  close(outfile)
  outfile=file(tempfilename,open="rt")
  text=readLines(outfile)
  close(outfile)
  unlink(tempfilename)
  if(mixtures==1)
     mixtures=NULL
  if(is.null(call))call=match.call()
  model = list(data = substitute(data), model = data$model,
        title = title, model.name = model.name, links = links, mixtures=mixtures,
        call = call, parameters=parameters,time.intervals=data$time.intervals, input = text, number.of.groups = number.of.groups,
        group.labels = group.labels, nocc = nocc, begin.time = data$begin.time, covariates=covariates,
        fixed=fixedvalues,design.matrix = complete.design.matrix, pims = pim,
        design.data = full.ddl,strata.labels=data$strata.labels,mlogit.list=mlogit.list)
  if(model.list$robust)model$nocc.secondary=nocc.secondary
#
#  If requested, simplify model which reconstructs PIMS, design matrix and rewrites the
#  MARK input file for the simplified model.
#
  model$profile.int=profile.int
  model$chat=chat
# Assign Mlogits that are set to a fixed value to a Logit link so they can be simplified
  if(mlogit0)
  {
	  fixedvalue=model$fixed$index
	  mlogit.indices=grep("mlogit",model$links)
	  if(length(mlogit.indices)>0 & length(fixedvalue)>0)
		  model$links[fixedvalue[fixedvalue%in%mlogit.indices]]="Logit"
  }
# Simplify the pim structure
  if(simplify) model=simplify.pim.structure(model)
#
# Check to make sure that the only rows in the design matrix that are all zeros are
# ones that correspond to fixed parameters.
#
  if(simplify)
  {                               
      dm=model$simplify$design.matrix
      fixed.rows=unique(model$simplify$pim.translation[model$fixed$index])
      zero.rows=(1:dim(dm)[1])[apply(dm,1,function(x) return(all(x=="0")))]
      if(length(fixed.rows)==0)
      {
         if(length(zero.rows)!=0)
            stop("One or more formulae are invalid because the design matrix has all zero rows for the following non-fixed parameters\n",
                  paste(row.names(dm)[zero.rows],collapse=","))
      }
      else
      {
         if(any(! (zero.rows%in%fixed.rows)))
            stop("One or more formulae are invalid because the design matrix has all zero rows for the following non-fixed parameters\n",
               paste(row.names(dm)[! (zero.rows%in%fixed.rows)],collapse=","))
      }
  }
#
#  Check to make sure that any parameter that used a sin link has an identity design matrix
#
  if(simplify)
  {
     for(i in 1:length(parameters))
     {
        parx=names(parameters)[i]
        if(model$parameters[[parx]]$link=="sin")
        {
           dm=model$simplify$design.matrix
           rows=unique(model$simplify$pim.translation[sort(unique(as.vector(unlist(sapply(model$pims[[parx]],function(x)x$pim[x$pim>0])))))])
           if(length(grep('[[:alpha:]]',as.vector(dm[rows,,drop=FALSE])))>0)
              stop("\nCannot use sin link with covariates")
           dm=suppressWarnings(matrix(as.numeric(dm),nrow=dim(dm)[1],ncol=dim(dm)[2]))
           if(any(rowSums(dm[rows,,drop=FALSE])>1) | any(colSums(dm[rows,,drop=FALSE])>1))
              stop("\nCannot use sin link with non-identity design matrix")
        }
     }
  }
#
#  check.mlogits(model)
#  model$mlogit.structure=NULL
  if(!is.null(model$simplify$links))
  {
     newlinks=model$simplify$links
     model$simplify$links=model$links
     model$links=newlinks
  }
  class(model)=c("mark",data$model)
  return(model)
}
