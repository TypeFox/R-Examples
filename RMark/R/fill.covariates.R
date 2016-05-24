#' Fill covariate entries in MARK design matrix with values
#' 
#' Replaces covariate names in design matrix with specific values to compute
#' estimates of real parameters at those values using the dataframe from
#' \code{\link{find.covariates}} after any value replacement.
#' 
#' The design matrix for a MARK model with individual covariates contains the
#' covariate names used in the model.  In computing the real parameters for the
#' encounter history of an individual it replaces instances of covariate names
#' with the individual covariate values.  This function replaces the cells in
#' the design matrix that contain individidual covariates with user-specified
#' values which is an edited version (if needed) of the dataframe returned by
#' \code{\link{find.covariates}}.
#' 
#' @param model MARK model object
#' @param values a dataframe matching structure of output from find.covariates
#' with the user-defined values entered
#' @return New design matrix with user-defined covariate values entered in
#' place of covariate names
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{find.covariates}}, \code{\link{compute.real}}
#' @keywords utility
#' @examples
#' \donttest{
#' data(dipper)
#' dipper$nsex=as.numeric(dipper$sex)-1
#' dipper$weight=rnorm(294)   
#' #NOTE:  This generates random valules for the weights so the answers using 
#' # ~weight will vary each time it is run
#' mod=mark(dipper,model.parameters=list(Phi=list(formula=~nsex+weight)))
#' # Show approach using individual calls to find.covariates, fill.covariates 
#' # and compute.real
#' fc=find.covariates(mod,dipper)
#' fc$value[fc$var=="nsex"]=0 # assign sex value to Female
#' design=fill.covariates(mod,fc) # fill design matrix with values
#' # compute and output survivals for females at average weight
#' female.survival=compute.real(mod,design=design)[1,] 
#' female.survival
#' # Next show same thing with a call to compute.real and a data frame for 
#' # females and then males
#' # compute and output survivals for females at average weight
#' female.survival=compute.real(mod,data=
#'       data.frame(nsex=0,weight=mean(dipper$weight)))[1,] 
#' female.survival
#' male.survival=compute.real(mod,data=data.frame(nsex=1,
#'          weight=mean(dipper$weight)))[1,] 
#' male.survival
#' # Fit model using sex as a group/factor variable and 
#' # compute v-c matrix for estimates
#' mod=mark(dipper,groups="sex",
#'      model.parameters=list(Phi=list(formula=~sex+weight)))
#' survival.by.sex=compute.real(mod,data=dipper,vcv=TRUE)
#' survival.by.sex$real[1:2]  # estimates
#' survival.by.sex$se.real[1:2] # std errors
#' survival.by.sex$vcv.real[1:2,1:2] # v-c matrix
#' survival.by.sex$vcv.real[1,2]/prod(survival.by.sex$se.real[1:2]) 
#' # sampling correlation of the estimates
#' }
fill.covariates <-
function(model,values)
{
# ----------------------------------------------------------------------------------------------------
#
#   fill.covariates  -   fills entries in design matrix with the values for covariates. The input argument 
#                        "values" is the output from find.covariates.
# Arguments:  
#
#   model    - MARK model object
#   values   - a dataframe matching structure of output from find.covariates
#
# Value:
#
#   xdesign  - new design matrix with user-defined covariate values entered.
#
# ----------------------------------------------------------------------------------------------------
model=load.model(model)
design=as.matrix(model$design.matrix)
values=values
if(!is.null(values))
for(i in 1:dim(values)[1])
   design[values$row[i],values$col[i]]=values$value[i]
xdesign=matrix(0,dim(design)[1],dim(design)[2])
for(i in 1:dim(design)[2])
   xdesign[,i]=as.numeric(design[,i])
return(xdesign)
}
