#' Merge time (occasion) and/or group specific covariates into design data
#' 
#' Adds new design data fields from a dataframe into a design data list
#' (\code{ddl}) by matching via time and/or group field in the design data.
#' 
#' Design data can be added to the parameter specific design dataframes with R
#' commands.  This function simplifies the process by enabling the merging of a
#' dataframe with a time and/or group field and one or more time and/or group
#' specific covariates into the design data list for a specific model
#' parameter. This is a replacement for the older function
#' \code{merge.occasion.data}. Unlike the older function, it uses the R
#' function \code{\link{merge}} but before merging it makes sure all of the
#' fields exist and that you are not merging data that already exists in the
#' design data.  It also maintains the row names in the case where design data
#' have been deleted prior to merging the design covariate data.
#' 
#' If \code{bytime=TRUE},the dataframe \code{df} must have a field named
#' \code{time} that matches 1-1 for each value of \code{time} in the design
#' data list (\code{ddl}).  All fields in \code{df} (other than time/group) are
#' added to the design data.  If you set \code{bygroup=TRUE} and have a field
#' named \code{group} in \code{df} and its values match the group fields in the
#' design data then group-specific values can be assigned for each time if
#' \code{bytime=TRUE}. If \code{bygroup=TRUE} and \code{bytime=FALSE} then it
#' matches by group and not by time.
#' 
#' @aliases merge_design.covariates merge.design.covariates
#' @param ddl current design dataframe for a specific parameter and not the
#' entire design data list (ddl)
#' @param df dataframe with time(occasion) and/or group-specific data
#' @param bygroup logical; if TRUE, then a field named \code{group} should be
#' in \code{df} and the values can then be group specific.
#' @param bytime logical; if TRUE, then a field named \code{time} should be in
#' \code{df} and the values can then be time specific.
#' @return Design dataframe (for a particular parameter) with new fields added.
#' See \code{\link{make.design.data}} for a description of the design data list
#' structure. The return value is only one element in the list rather than the
#' entire list as with the older function \code{merge.occasion.data}.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{make.design.data}}, \code{\link{process.data}}
#' @keywords utility
#' @examples
#' 
#' data(dipper)
#' dipper.proc=process.data(dipper)
#' ddl=make.design.data(dipper.proc)
#' df=data.frame(time=c(1:7),effort=c(10,5,2,8,1,2,3))
#' # note that the value for time 1 is superfluous for CJS but not for POPAN
#' # the value 10 will not appear in the summary because there is no p for time 1
#' summary(ddl$p)
#' ddl$p=merge_design.covariates(ddl$p,df)
#' summary(ddl$p)
#' #Statement below will create an error because a value for time 7 not given
#' #ddl=merge.occasion.data(dipper.proc,ddl,"p",data.frame(time=c(1:6),effort=c(10,5,2,8,1,2)))
#' #
#' # Assign group-specific values
#' #
#' 
#' data(dipper)
#' dipper.proc=process.data(dipper)
#' ddl=make.design.data(dipper.proc)
#' df=data.frame(time=c(1:7),effort=c(10,5,2,8,1,2,3))
#' # note that the value for time 1 is superfluous for CJS but not for POPAN
#' # the value 10 will not appear in the summary because there is no p for time 1
#' summary(ddl$p)
#' ddl$p=merge_design.covariates(ddl$p,df)
#' summary(ddl$p)
#' #Statement below will create an error because a value for time 7 not given
#' #ddl=merge.occasion.data(dipper.proc,ddl,"p",data.frame(time=c(1:6),effort=c(10,5,2,8,1,2)))
#' #
#' # Assign group-specific values
#' #
#' dipper.proc=process.data(dipper,groups="sex")
#' ddl=make.design.data(dipper.proc)
#' df=data.frame(group=c(rep("Female",6),rep("Male",6)),time=rep(c(2:7),2),
#' 		effort=c(10,5,2,8,1,2,3,20,10,4,16,2))
#' merge_design.covariates(ddl$p,df,bygroup=TRUE)
#' 
"merge_design.covariates"<-function(ddl,df,bygroup=FALSE,bytime=TRUE)
#                   - enables design data fields to be added to the design.data(ddl)
#                    that match time fields in design data (eg effort), by group (eg site covariates)
#                    or by time and group - site and time dependent covariates.
#
# Arguments:
#
# ddl              - current design dataframe
# df               - dataframe with time and/or group specific data
# bygroup          - if TRUE, values are group specific
# bytime           - if TRUE, values are time-specific
#
# Value:
#
#  ddl - modified design data
#
# -------------------------------------------------------------------------------------------------------------
{
#
# Check to make sure bytime or bygroup
#
  if(!bygroup & !bytime) stop("\n either bygroup or bytime must be TRUE\n")
#
# Check to make sure df is a dataframe, contains the field time (if bytime=TRUE) and
# that a value is given for each time. Likewise do the same for group.
#
  if(!is.data.frame(df))stop(paste("\n",substitute(df),"is not a dataframe"))
  if(bytime)
  {
     if(is.null(df$time))stop(paste("\n",substitute(df),"does not contain field named time"))
     if(any(!ddl$time%in%unique(df$time)))
         stop(paste("\n",substitute(df),"does not contain a time value for each time in design data"))
     if(!bygroup & !is.null(df$group)) stop(paste("\n",substitute(df),"contains a group field but bygroup=FALSE"))
  }
  if(bygroup)
  {
     if(is.null(df$group))stop(paste("\n",substitute(df),"does not contain field named group"))
     if(any(!ddl$group%in%unique(df$group)))
         stop(paste("\n",substitute(df),"does not contain a group value for each group in design data"))
     if(!bytime & !is.null(df$time)) stop(paste("\n",substitute(df),"contains a time field but bytime=FALSE"))
  }
  if(bygroup & bytime)
  {
     all.gt=paste(ddl$group,ddl$time,sep="")
     all.gt.df=paste(df$group,df$time,sep="")
     if(any(!all.gt%in%all.gt.df))
         stop(paste("\n",substitute(df),"does not contain a group/time value for each group/time in design data"))
  }
#
# Check to make sure there is no overlap in field names other than group & time
#
  if(any(names(ddl)%in% names(df)[!names(df)%in%c("time","group")]))
      stop(paste("\n",substitute(df),"uses the same field names as used in design data\n"))
#
# Save row names for the case in which design data have been deleted
#
   save.row.names=row.names(ddl)
#
# Add sequence field to resort after merge
#
  ddl$xxsequence=1:dim(ddl)[1]
#
# Merge data
#
  ddl=merge(ddl,df,sort=FALSE)
#
# Re-sort data and remove field
#
  ddl=ddl[order(ddl$xxsequence),]
  ddl$xxsequence=NULL
#
# reassign row names and return data
#
  row.names(ddl)=save.row.names
  return(ddl)
}
