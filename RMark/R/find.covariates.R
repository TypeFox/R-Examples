#' Find covariates in MARK design matrix
#' 
#' Finds and extracts cells in MARK design matrix containing covariates.
#' Computes mean values of the covariates and assigns those as default values.
#' Returns dataframe that can be edited to replace default values which are
#' then inserted into the design matrix with \code{\link{fill.covariates}} to
#' enable computation of estimates of real parameters with
#' \code{\link{compute.real}}.
#' 
#' The design matrix for a MARK model with individual covariates contains
#' entries with the covariate names used in the model. In computing the real
#' parameters for the encounter history of an individual it replaces instances
#' of covariate names with the individual covariate values.  This function
#' finds all of the cells in the design matrix that contain individidual
#' covariates and constructs a dataframe of the name of the real parameter, the
#' position (row, col) in the design matrix and a default value for the
#' covariate. The default field value is assigned to one of three values in the
#' following priority order: 1) the mean value for the covariates in data (if
#' data is not NULL), 2) the mean values used in the MARK output (if
#' data=NULL,usemean=TRUE), 3) 0 (if usemean=FALSE and data=NULL). The values
#' can also be modified using \code{fc=edit(fc)} where \code{fc} is the value
#' from this function.
#' 
#' @param model MARK model object
#' @param data dataframe used to construct MARK model object; not processed
#' data list
#' @param usemean logical; if TRUE uses mean value of covariate for default and
#' otherwise uses 0
#' @return A dataframe with the following fields \item{rnames}{name of real
#' parameter} \item{row}{row number in design matrix (equivalent to
#' \code{parm.indices} in call to \code{\link{compute.real}}} \item{col}{column
#' number in design matrix} \item{var}{name of covariate} \item{value}{value
#' for covariate}
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{fill.covariates}}, \code{\link{compute.real}}
#' @keywords utility
#' @examples
#' 
#' # see examples in fill.covariates
#' 
find.covariates <-
function(model,data=NULL,usemean=TRUE)
{
# ----------------------------------------------------------------------------------------------------
#
#   find.covariates  - finds entries in design matrix that contain covariate values.  It creates a dataframe
#                      with the rownames, row and col #s and the covariate name.  The field value is assigned
#                      to one of three values in the following priority order
#                         1) the mean value for the covariates in data (if data is not NULL)
#                         2) the mean values used in the MARK output (if data=NULL,usemean=TRUE)
#                         3) 0 (if usemean=FALSE and data=NULL).
#                      If data is specified but some covariates are skipped it will use #2 or #3 to
#                      fill in the remainder of the covariates.
#                      The values output here can also be modified using
#                      fc=edit(fc) where fc is the output from this function, but
#                      use of data is better for scripting (automating) the analysis.
#
# Arguments:  
#   model   - MARK model object
#   data    - dataframe of variables which are averaged for prediction
#   usemean - if TRUE, assigns mean of covariate value as default value; otherwise uses 0 
#
# Value:
#   cov.datalist  - a dataframe with fields: row, col, var, and value.
#
# ----------------------------------------------------------------------------------------------------
model=load.model(model)
design=model$design.matrix
#
#  If usemean=TRUE compute means of numeric covariates; otherwise set to zero
#
if(!is.null(model$covariates))
{
   if(!usemean)
   {
      means=rep(0,length(model$covariates))
   }
   else
   {
     if(is.null(model$results$covariate.values))
     {
        means=rep(0,length(model$covariates))
        warning("No value given for data argument so any covariate values are set to 0.")
     }
     else
     {
        means=model$results$covariate.values$Value
     }
   }
   names(means)=model$covariates
}
#
#  This will only be reached if this function is used on old model results that
#  don't have a model$covariates element
#
else
   if(is.null(data))
      stop("\n This is an old model object. Re-run model or use the data argument must be specified and be a dataframe. Use data and not processed data list.\n")
   else
#
#  If there are no covariates in the model return NULL
#
      return(NULL)
if(!is.null(data))
{
   if(!is.data.frame(data))
      stop("\n data argument must be a dataframe. Use data and not processed data list.\n")
   data.means=apply(as.matrix(data[,sapply(data,is.numeric),drop=FALSE]),2,mean)
   coln=names(means)%in%names(data.means)
#   coln=match(names(data.means),names(means))
   colnd=match(names(means),names(data.means))
   means[coln]=data.means[colnd[!is.na(colnd)]]
#   means[coln[!is.na(coln)]]=data.means[colnd[!is.na(colnd)]]
}
#
# Next create a temp environment and assign variables with covariate values
# and create the product function to compute product(x,y) entries in the
# design matrix
#
temp.env=new.env()
for (i in 1:length(means))
  assign(names(means)[i],means[i],envir=temp.env)
product=function(x,y)return(x*y)
assign("product",product,envir=temp.env)
#
# Find the columns with covariates (covcol)
#
values=NULL
cov.datalist=NULL
cov.rownames=NULL
covcol=apply(design,2,function(x){any(is.na(suppressWarnings(as.numeric(x))))})
covcol=(1:length(covcol))[covcol]
#
# Loop over each column with a covariate and examine each row entry.  If it is
# not numeric evaluate it's value in the temporary environment and store the
# value, row and column number and its name.
#
for(j in covcol)
{
   for(i in 1:dim(design)[1])
   {
      if(is.na(suppressWarnings(as.numeric(design[i,j]))))
      {
         cov.rownames=c(cov.rownames,rownames(design)[i])
         cov.datalist=rbind(cov.datalist,c(i,j,design[i,j]))
         values=c(values,eval(parse(text=design[i, j]),envir=temp.env))
      }
   }
}
#
# If there were no columns with covariates, return NULL otherwise return
# a dataframe with the position and covariate value which can be used in fill.covariates
# to fill the design matrix with the assigned value.
#
   if(is.null(cov.rownames))
      return(NULL)
   else
      return(cov.datalist=data.frame(rnames=cov.rownames,row=as.numeric(cov.datalist[,1]),
             col=as.numeric(cov.datalist[,2]),var=cov.datalist[,3],value=values))
}
