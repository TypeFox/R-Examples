#' Make time-varying dummy variables from time-varying factor variable
#' 
#' Create a new dataframe with time-varying dummy variables from a time-varying
#' factor variable.  The time-varying dummy variables are named appropriately
#' to be used as a set of time dependent individual covariates in a parameter
#' specification
#' 
#' An example of the var.name and times is var.name="observer", times=1:5. The
#' code expects to find observer1,...,observer5 to be factor variables in x. If
#' there are k unique levels (excluding ".") across the time varying factor
#' variables, then k-1 dummy variables are created for each of the named factor
#' variables.  They are named with var.name, level[i], times[j] concatenated
#' together where level[i] is the name of the facto level i.  If there a m
#' times then the new data set will contain m*(k-1) dummy variables. If the
#' factor variable includes any "." values these are ignored because they are
#' used to indicate a missing value that is paired with a missing value in the
#' encounter history. Note that it will create each dummy variable for each
#' factor even if a particular level is not contained within a factor (eg
#' observers 1 to 3 used but only 1 and 2 on occasion 1).
#' 
#' @param x dataframe containing set of factor variables with names composed of
#' var.name prefix and times suffix
#' @param var.name prefix for variable names
#' @param times numeric suffixes for variable names
#' @param intercept the value of the factor variable that will be used for the
#' intercept
#' @param delete if TRUE, the origninal time-varying factor variables are
#' removed from the returned dataframe
#' @return x: a dataframe containing the original data (with time-varying
#' factor variables removed if delete=TRUE) and the time-varying dummy
#' variables added.
#' @author Jeff Laake
#' @export
#' @keywords utility
#' @examples
#' 
#' # see example in weta
#' 
make.time.factor=function(x,var.name,times,intercept=NULL,delete=TRUE)
#
#  This function takes a time varying factor variable and creates a
#  a set of time varying dummy variables to be used as time dependent
#  individual covariates.
#
#  Arguments:
#      x        - dataframe
#     var.name  - prefix variable name
#     times     - suffix numbers for variable names
#     intercept - value to used for intercept
#     delete    - if TRUE remove the original columns from the data
#
#  Value:
#     dataframe with new time varying dummy variables added
#
{
#  create var.names from prefix and times and check to make sure they exist
   var.names=paste(var.name,times,sep="")
   if(any(!var.names %in%names(x))) stop(paste("Following fields not found: ",
           paste(var.names[!var.names %in%names(x)],collapse=",")))
#  create initial new dataframe and delete orginal columns if requested
   if(delete)
      y=subset(x,select=names(x)[!names(x)%in%var.names])
   else
      y=x
# Compute unique levels across all of the factor variables
  ulevels=NULL
  for(i in 1:length(var.names))ulevels=c(ulevels,levels(x[,var.names[i]]))
  ulevels=unique(ulevels)
# Loop over each time-varying factor variable and create k-1 time varying dummy
# variables for each factor variable where k is the number of levels of the factor
# which is the length of times.  Any "." values are ignored.
#
   for(i in 1:length(var.names))
   {
     xx=x[,var.names[i],drop=FALSE]
     char=unlist(strsplit(var.names[i],""))
     startc=(1:length(char))*as.numeric(char%in%0:9)
     startc=min(startc[startc!=0])
     time.index=substr(var.names[i],start=startc,stop=length(char))
     names(xx)="x"
     xx$x=factor(as.character(xx$x),levels=ulevels)
     mat=data.frame(model.matrix(~-1+x,xx))
     mat=subset(mat,select=names(mat)[names(mat)!="x."])
     xlevels=ulevels
     if(!is.null(intercept))
     {
        mat=subset(mat,select=names(mat)[names(mat)!=paste("x",intercept,sep="")])
        xlevels=ulevels[ulevels!=intercept]
     }
     names(mat)=paste(substr(var.name,start=1,stop=startc-1),xlevels[xlevels!="."],time.index,sep="")
     y=cbind(y,mat)
   }
   return(y)
}
