#' Export output files for appending into MARK .dbf/.fpt format
#' 
#' Creates renamed versions of the output,vcv and residual files so they can be
#' appended into a MARK .dbf file.
#' 
#' If \code{model} is a marklist then it exports each model in the marklist.
#' The function simply copies the files with new names so the MARK interface
#' will recognize them. The marknnn.out is copied as marknnnY.tmp, marknnn.res
#' is copied as marknnnx.tmp and marknnn.vcv is copied as marknnnV.tmp.  You
#' can create a MARK .dbf by using \code{\link{export.chdata}} to create an
#' input file for MARK, opening MARK (MARKINT.EXE) to create a new .dbf with
#' the input file, and then using the Output/Append to select the output file
#' (marknnnY.tmp) to append the model with its files. Then you can use any
#' facilities of MARK that are not already included in RMark.
#' 
#' ***Warning*** Make sure that you use the .inp created by
#' \code{\link{export.chdata}} with your processed data to create the MARK .dbf
#' file rather than using a separate similar .inp file.  It is essential that
#' the group structure and ordering of groups matches between the .inp file and
#' the exported models or you can get erroneous results.
#' 
#' @param model a mark model object or marklist object
#' @param replace if file exists and replace=TRUE, file will be over-written
#' @return None
#' @author Jeff Laake
#' @seealso \code{\link{export.chdata}}
#' @keywords utility
#' @export
#' @examples
#' 
#' data(dipper)
#' mymodel=mark(dipper,threads=1)
#' export.model(mymodel,replace=TRUE)
#' 
export.model <-
function(model,replace=FALSE)
{
# -----------------------------------------------------------------------------------------------------------------------
#
# export.model   -   creates files for appending into MARK .dbf/.fpt files
#
# Arguments:
#
# model            - a mark model object or marklist of models
# replace          - if file exists and replace=TRUE, file will be over-written
#
# Value:
#
#  none
#
#
  model=load.model(model)
  if(class(model)[1]=="marklist")
     for(i in 1:(length(model)-1))
        export.model(model[[i]],replace=replace)
  else
  {
     file.copy(paste(model$output,".vcv",sep=""),paste(model$output,"V.tmp",sep=""), overwrite = replace)
     file.copy(paste(model$output,".out",sep=""),paste(model$output,"Y.tmp",sep=""), overwrite = replace)
     file.copy(paste(model$output,".res",sep=""),paste(model$output,"X.tmp",sep=""), overwrite = replace)
  }
  invisible()
}
