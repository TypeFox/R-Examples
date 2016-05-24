#' Export capture-history data to MARK .inp format
#' 
#' Creates a MARK .inp file from processed data list that can be used to create
#' a MARK .dbf file for use with MARK directly rather than with the RMark
#' package.
#' 
#' The default is to include none of the covariates in the processed data list.
#' All of the covariates can be passed by setting covariates="all".
#' 
#' After you have created the MARK .dbf/.fpt files with the exported .inp file,
#' then you can use \code{\link{export.chdata}} to export models that can be
#' imported into the MARK interface.  However note the following: ***Warning***
#' Make sure that you use the .inp created by \code{\link{export.chdata}} with
#' your processed data to create the MARK .dbf file rather than using a
#' separate similar .inp file.  It is essential that the group structure and
#' ordering of groups matches between the .inp file and the exported models or
#' you can get erroneous results.
#' 
#' @param data processed data list resulting from process.data
#' @param filename quoted filename (without .inp extension)
#' @param covariates vector of names of covariate variables in data to include
#' @param replace if file exists and replace=TRUE, file will be over-written
#' @return None
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{import.chdata}}
#' @keywords utility
#' @examples
#' 
#' data(dipper)
#' dipper$numeric.sex=as.numeric(dipper$sex)-1
#' dipper.processed=process.data(dipper,group="sex")
#' export.chdata(dipper.processed, filename="dipper", 
#'          covariates="numeric.sex",replace=TRUE)
#' #
#' # Had sex been used in place of numeric.sex in the above command, 
#' # MARK would have been unable to use it as a covariate
#' # because it is not a numeric field
#' 
export.chdata <-
function(data, filename, covariates=NULL, replace=FALSE)
{
# -----------------------------------------------------------------------------------------------------------------------
#
# export.chdata   -   creates a MARK .inp file from dataframe that can be used to create a MARK .dbf file
#
# Arguments:
#
# data             - processed data list resulting from process.data
# filename         - filename (without .inp extension)
# covariates       - names of covariate variables in data to put in .inp; default is to use none and
#                    all can be passed by setting covariates="all")
# replace          - if file exists and replace=TRUE, file will be over-written
#
# Value:
#
#  none; creates text file 
#
#
#
  outfile=paste(filename,".inp",sep="")
  if(file.exists(outfile))
     if(replace)
     {
        file.create(outfile)
     }
     else
        stop("File already exists and replace=FALSE")
#
# Check to make sure this is a processed data list
#
  if(is.null(data$data))
     stop("\nUse processed data list and not original dataframe for the data argument\n")
#
# Output data portion of MARK input file:
#
  if(data$model!="Nest")
  {
     ch=data$data$ch
     zz=as.data.frame(ch)
     zz=cbind(zz,data$freq)
     if(!is.null(covariates))
     {
        if(covariates[1]=="all")
           zz=data.frame(cbind(zz,data$data[,-1])) 
        else
           zz=data.frame(cbind(zz,data$data[,covariates])) 
     }
#  
#    This outputs capture history, frequency and any covariates
#
	 write.table(zz,file=outfile,eol=";\n",sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
 } else
#	 
#    Output nest survival model
#
 {
	 
	 zz=data$data[,1:5]
	 if(!is.null(covariates))
	 {
		 if(covariates[1]=="all")
			 zz=data.frame(cbind(zz,data$data[,-(1:5)])) 
		 else
			 zz=data.frame(cbind(zz,data$data[,covariates])) 
	 }
	 if(is.null(data$group.covariates)) ng=1 else ng=nrow(data$group.covariates)
	 for(i in 1:ng)
	 {
		write(paste("Nest survival group =",i,";"),file=outfile,append=TRUE)
		write.table(zz[zz$group==i,],file=outfile,eol=";\n",sep=" ",col.names=FALSE,row.names=FALSE,quote=FALSE,append=TRUE)
	 }	 
  }
  invisible()
}
