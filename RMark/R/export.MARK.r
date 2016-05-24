#' Export data and models for import in MARK
#' 
#' Creates a .Rinp, .inp and optionally renamed output files that can be
#' imported into MARK to create a new MARK project with all the data and output
#' files.
#' 
#' If you use Nest model and NestAge covariate you should use the processed
#' data list (model$data) from a model using NestAge in the formula because the
#' necessary individual covariates are added to the processed data list.  Also,
#' use default of ind.covariates="all".
#' 
#' After running this function to export the data and models to the working
#' directory, start program MARK and select the File/RMARK Import menu item.
#' Navigate to the working directory and select the "project.name".Rinp file.
#' MARK will take over and create the project files and will import any
#' specified model output files.
#' 
#' DO NOT use a "project name" that is the same as a filename you use as input
#' (*.inp) to RMark because this function will create a file "project name".inp
#' and it would over-write your existing file and the format could change. If
#' you try this, the code will return an error that the filename is invalid
#' because the file already exists.
#' 
#' @param x processed data list used to build models
#' @param project.name character string to be used for prefix of filenames and
#' for MARK project name; do not use "." in the filename
#' @param model either a single mark model or a marklist
#' @param replace if TRUE it will replace any existing files
#' @param chat user-specified chat value if desired
#' @param title MARK project title string
#' @param ind.covariates vector of character strings specifying names of
#' individual covariates or "all" to use all present
#' @return None
#' @export
#' @author Jeff Laake
#' @seealso \code{\link{export.chdata}}, \code{\link{export.model}}
#' @keywords utility
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' data(mallard)
#'Dot=mark(mallard,nocc=90,model="Nest",
#'		model.parameters=list(S=list(formula=~1)))
#'mallard.proc=process.data(mallard,nocc=90,model="Nest")
#'export.MARK(mallard.proc,"mallard",Dot,replace=TRUE)
#'data(robust)
#'time.intervals=c(0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0)
#'S.time=list(formula=~time)
#'p.time.session=list(formula=~-1+session:time,share=TRUE)
#'GammaDoublePrime.random=list(formula=~time,share=TRUE)
#'model.1=mark(data = robust, model = "Robust",
#'		time.intervals=time.intervals,
#'		model.parameters=list(S=S.time,
#'				GammaDoublePrime=GammaDoublePrime.random,p=p.time.session))
#'robust.proc=process.data(data = robust, model = "Robust",
#'		time.intervals=time.intervals)
#'export.MARK(robust.proc,"robust",	model.1,replace=TRUE)
#'}
export.MARK=function(x,project.name,model=NULL,replace=FALSE,chat=1.0,title="",ind.covariates="all")
{
# exports model and data files for import into the MARK interface
#
# Arguments:
# 
# x              - processed dataframe used to build models
# project.name   - character string to be used for prefix of filenames and for MARK project name
# model          - either a single mark model or a marklist
# replace        - if TRUE it will replace any existing files
# chat           - user-specified chat value if desired
# title          - MARK project title string
# ind.covariates - vector of character strings specifying names of individual covariates
#
# Value: no returned value, but it will create a "project.name".inp,"project.name".Rinp, 
#        and copy the output files for the mark models to the names needed to import into MARK
#
  if(is.null(x$model)) stop(paste(substitute(x), "is not a processed dataframe"))
  filename=paste(project.name,".Rinp",sep="")
  if(file.exists(filename)&!replace) 
  {
    stop("Project already exists and replace=FALSE")
  } else
  {
    if(!file.exists(filename))
    { 
       xfilename=paste(project.name,".inp",sep="")
       if(file.exists(xfilename)) stop("Use a different project name because the .inp file already exists")
    }
  }
  if(x$model=="MultScalOcc")
  {
	 x$mixtures= match("1",x$time.intervals)
	 x$time.intervals=rep(1,length(x$time.intervals))
	 x$nocc=nchar(x$data$ch[1])
  }
  if(x$model=="Nest") 
	  nocc=1
  else
	  nocc=nchar(x$data$ch[1])
  write(setup.model(x$model, nocc, x$mixtures)$etype,file=filename)
  write(x$mixtures,file=filename,append=TRUE)
  if(setup.model(x$model, nocc, x$mixtures)$robust)
     write(nocc,file=filename,append=TRUE)
  else
	  write(x$nocc,file=filename,append=TRUE)
  if(is.null(ind.covariates))
  {
    write("0",file=filename,append=TRUE)  
  } else
  {
	if(ind.covariates[1]=="all")
	{
		ind.covariates=names(x$data[!sapply(x$data,is.factor)])
		if(x$model=="Nest")
			ind.covariates=ind.covariates[!ind.covariates%in%c("FirstFound","LastPresent","LastChecked","Fate","freq")]
		else
			ind.covariates=ind.covariates[!ind.covariates%in%c("ch","freq")]	
	}
	if(length(ind.covariates)==0)
	{
		ind.covariates=NULL
		write("0",file=filename,append=TRUE)  
	} else
	{
       if(all(ind.covariates%in%names(x$data)))
       {
          write(length(ind.covariates),file=filename,append=TRUE)
          write(matrix(ind.covariates,ncol=1),file=filename,append=TRUE)
       } else
          stop("One or more of the ind.covariates values were not found in the data")
    } 
  }
  write(x$nstrata,file=filename,append=TRUE)
  if(is.null(x$strata.labels))
    x$strata.labels=rep("NA",x$nstrata)
  write(matrix(paste(x$strata.labels,":Stratum",1:length(x$strata.labels),sep=""),ncol=1),file=filename,append=TRUE)
  write(x$time.intervals,ncolumns=length(x$time.intervals),file=filename,append=TRUE)
  if(is.null(x$group.covariates))
    write("1",file=filename,append=TRUE)  
  else    
    write(nrow(x$group.covariates),file=filename,append=TRUE)
  write(matrix(colnames(x$freq),ncol=1),file=filename,append=TRUE)
  write(nrow(x$freq),file=filename,append=TRUE)
  write(chat,file=filename,append=TRUE)
  if(title=="")
    write("NA",file=filename,append=TRUE)
  else
    write(title,file=filename,append=TRUE)
  export.chdata(x,project.name,covariates=ind.covariates,replace=replace)
  if(!is.null(model)) export.model(model,replace=replace)
  return(NULL)
}



