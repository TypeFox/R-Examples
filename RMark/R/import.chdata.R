#' Import capture-recapture data sets from space or tab-delimited files
#' 
#' A relatively flexible function to import capture history data sets that
#' include a capture (encounter) history read in as a character string and an
#' arbitrary number of user specified covariates for the analysis.
#' 
#' This function was written both to be a useful tool to import data and as an
#' example for more specific import functions that a user may want to write for
#' data files that do not satisfy the requirements of this function.  In
#' particular this function will not handle files with fixed-width format files
#' that do not contain appropriate tab or space delimiters between the fields.
#' It also requires that the first field is the capture (encounter) history
#' which is named "ch" and is a character string.  The remaining fields are
#' arbitrary in number and type and are user defined based on the arguments to
#' the functions. Variables that will be used for grouping should be defined
#' with the \code{field.type="f"}. Numeric individual covariates (e.g., weight)
#' should be input as \code{field.type="n"}. Fields in the file that should not
#' be imported should be assigned \code{field.type="s"}.  The examples below
#' illustrate different uses of the calling arguments to import several
#' different data sets that meet the modest requirements of this function.
#' 
#' If you specify a frequency for the encounter history, the field name must be
#' \code{freq}.  If you use any other name or spelling it will not be
#' recognized and the default frequency of 1 will be used for each encounter
#' history.  This function should not be used with files structured for input
#' into the MARK interface.  To use those types of files, see
#' \code{\link{convert.inp}}.  It is not neccessary to use either function to
#' create a dataframe for RMark.  All you need to is create a dataframe that
#' meets the specification of the RMark format.  For example, if you are
#' simulating data, you only need to create a dataframe with the fields ch,
#' freq (if differs from 1) and any covariates you want and then you can use
#' \code{\link{process.data}} on the dataframe.
#' 
#' If you have comments in your data file, they should not have a column header
#' (field name in first row).  If \code{use.comments=TRUE} the comments are
#' used as row names of the data frame and they must be unique.  If
#' \code{use.comments=FALSE} and the file contains comments they are stripped
#' out.
#' 
#' @param filename file name and path for file to be imported; fields in file
#' should be space or tab-delimited
#' @param header TRUE/FALSE; if TRUE first line is name of variables
#' @param field.names vector of field names if header=FALSE; first field should
#' always be ch - capture history remaining number of fields and their names
#' are arbitrary
#' @param field.types vector identifying whether fields (beyond ch) are numeric
#' ("n") or factor ("f") or should be skipped ("s")
#' @param use.comments if TRUE values within /* and */ on data lines are used
#' as row.names for the RMark dataframe.  Only use this option if they are
#' unique values.
#' @return A dataframe for use in MARK analysis with obligate \code{ch}
#' character field representing the capture (encounter) history and optional
#' covariate/grouping variables.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{export.chdata}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' pathtodata=paste(path.package("RMark"),"extdata",sep="/")
#' example.data<-import.chdata(paste(pathtodata,"example.data.txt",sep="/"),
#'       field.types=c("n","f","f","f"))
#' edwards.eberhardt<-import.chdata(paste(pathtodata,"edwardsandeberhardt.txt",
#'       sep="/"),field.names="ch",header=FALSE)
#' dipper<-import.chdata(paste(pathtodata,"dipper.txt",sep="/"),
#'       field.names=c("ch","sex"),header=FALSE)
#' }
import.chdata <-
function(filename, header=TRUE, field.names=NULL, field.types=NULL, use.comments=TRUE)
{
#
# import.chdata - reads in capture history and user specified covariates for analysis of 
#                 mark-recapture data. 
#
# Arguments:
#
# filename      - file name and path for file to be imported; fields in file should be space or tab-delimited
# header        - TRUE/FALSE; if TRUE first line is name of variables
# fied.names    - vector of field names if header=FALSE; first field should always be ch - capture history
#                 remaining number of fields and their names are arbitrary 
# field.types   - vector identifying whether fields (beyond ch) are numeric ("n") or factor ("f") or should be skipped ("s")
# use.comments  - logical; if TRUE values within /* and */ on data lines are
#                            used as row.names for the RMark dataframe.  Only use if
#                            they are unique values.
#
#
# Value: dataframe for use in MARK analysis with obligate ch field and optional covariate/grouping variables
#
#
strip.list=strip.comments(filename,use.comments=use.comments,header=header)
rn=strip.list$rn
out.filename=strip.list$out.filename
filename=out.filename
if(!is.null(field.names))header=FALSE
data=read.table(filename,colClasses=c("character"),header=header)
unlink(out.filename)
if(header & names(data)[1]!="ch") stop("First field should be named ch; Either first row doesn't contain field names or first field not named properly")
nvar=dim(data)[2]
#
# Assign field names if they are not the first line of the data file
#
if(!header)
  if(nvar==length(field.names))
     if(field.names[1]=="ch")
        names(data)=field.names
     else
        stop("First field should be the capture-history and named ch in field.names")
  else
     stop("Length of field.names does not match number of columns in data")
#
# If field.types are specified create factor and numeric variables as assigned
# Otherwise presume they are all factors
#
if(nvar>1)
{
   if(is.null(field.types)) field.types=rep("f",nvar-1)
   for( i in 2:dim(data)[2])
   {
       if(field.types[i-1] =="f")
          data[,i]=as.factor(data[,i])
       else
          if(field.types[i-1]=="n")
             data[,i]=as.numeric(data[,i]) 
   }
   data=data[,!c(FALSE,field.types=="s")]
}
row.names(data)=rn
return(data)
}
