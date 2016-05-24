#' Removes unused MARK output files
#' 
#' Identifies all unused (orphaned) mark*.inp, .vcv, .res and .out and .tmp
#' files in the working directory and removes them.  The orphaned files are
#' determined by examining all mark objects and lists of mark objects (created
#' by \code{\link{collect.models}}) to create a list of files in use. All other
#' files are treated as orphans to delete.
#' 
#' This function removes orphaned output files from MARK.  This occurs when
#' there are output files in the subdirectory that are not associated with a
#' mark object in the current R session (.Rdata). For example, if you repeat an
#' analysis or set of analyses and store them in the same object then the
#' original set of output files would no longer be linked to an R object and
#' would be orphaned.
#' 
#' As an example, consider the \code{\link{mallard}} analysis. The first time
#' you run the analysis script in an empty subdirectory it would create 9 sets
#' of MARK output files (mark001.out,.vcv,.res,.inp to
#' mark009.out,.vcv,.res,.inp) and each would be linked to one of the objects
#' in \code{nest.results}. When the command
#' \code{AgePpnGrass=nest.results$AgePpnGrass} was issued, both of those
#' \code{mark} objects were linked to the same set of output files.  Now if you
#' were to repeat the above commands and re-run the models and stored the
#' results into \code{nest.results} again, it would create files with prefixes
#' 10 to 18. Because that would have replaced \code{nest.results}, none of the
#' files numbered 1 to 9 would be linked to an object.
#' \code{cleanup(ask=FALSE)} automatically removes those orphan files.  If you
#' delete all objects in the R session with the command
#' \code{rm(list=ls(all=TRUE))}, then subsequently \code{cleanup(ask=FALSE)}
#' will delete all MARK output files because all of them will be orphans.
#' Output files can also become orphans if MARK finishes but for some reason R
#' crashes or you forget to save your session before you exit R.  Orphan output
#' files can be re-linked to an R object without re-running MARK by re-running
#' the \code{\link{mark}} function in R and specifing the \code{filename}
#' argument to match the base portion of the orphaned output file (eg
#' "mark067").  It will create all of the necessary R objects and then asks if
#' you want to use the existing file.  If you respond affirmatively it will
#' link to the orphaned files.
#' 
#' @param lx listing of R objects; defaults to list of workspace from calling
#' environment; if NULL it uses ls(envir=parent.frame())
#' @param ask if TRUE, prompt whether each file should be removed. Typically
#' will be used with ask=FALSE but default of TRUE may avoid problems
#' @param prefix prefix for filename if different than "mark"
#' @return None
#' @export
#' @author Jeff Laake
#' @keywords utility
cleanup <- function(lx=NULL,ask=TRUE,prefix="mark")
# ----------------------------------------------------------------------------------------
#
# cleanup   - remove unused MARK* files
#
# Arguments:
#
#  lx             - listing of R objects
#  ask            - if TRUE, ask whether each file should be deleted
#
# Value:
#
#  None
#
#
#  10 Jan 06; added to cleanup files left by deleted mark objects
# ----------------------------------------------------------------------------------------
{
purge=function(type,model.filenames,ask)
{
   file.list=list.files(pattern=paste(prefix,"[0123456789]+[:.:]",type,sep=""))
   if(length(file.list)>0)
   for (i in 1:length(file.list))
   {
     if(!file.list[i] %in% paste(model.filenames,".",type,sep=""))
     {
        if(ask)
           answer=readline(paste("Delete file",file.list[i],"(y/n)[y]?"))
        else
           answer="y"
        if(substr(answer,1,1)!="n")
           unlink(file.list[i])
     }
   }
}
#
# Collect mark model objects
#
xx="\\\\."
if(is.null(lx))lx=ls(envir=parent.frame())
exclude=grep("\\*",lx)
if(length(exclude)>0)lx=lx[-exclude]
model.list=collect.model.names(lx,warning=FALSE)
model.filenames=NULL
if(!is.null(model.list))
{
   for( i in 1:length(model.list))
   {
      if(eval(parse(text=paste("is.list(",model.list[i],")",sep=""))))
         model.filenames=c(model.filenames,eval(parse(text=paste(model.list[i],"$output",sep=""))))
      else
      {
         zz=eval(parse(text=model.list[i]))
         model.filenames=c(model.filenames,eval(parse(text=paste('strsplit("',zz,'",','"',xx,'"',")[[1]][1]",sep=""))))
      }
   }
}        
#  x in the function call did not work in Linux; used xzx and worked fine
#myf=function(xzx)
#{
#   eval(parse(text=paste("ifelse(is.list(xzx),xzx$output,strsplit(xzx,'",xx,"')[[1]][1])",sep="")))
#}
#debug(myf)
#
#model.filenames=sapply(model.list,myf)
#model.filenames=sapply(paste(model.list,"$output",sep=""),function(xzx)eval(parse(text=xzx)))
#
# Collect marklist objects
#
model.list=NULL
for (i in 1:length(lx))
{
   classval=class(eval(parse(text=lx[i]),envir=parent.frame(2)))
   if(classval[1]=="marklist")
   model.list=c(model.list,lx[i])
}
blank=""
if(length(model.list)!=0)
   for(i in 1:length(model.list))
   {
      ml=eval(parse(text=model.list[i]),envir=parent.frame(2))
      num.models=length(ml)-as.numeric(!is.null(ml$model.table))
      for (j in 1:num.models)
         model.filenames=c(model.filenames,
           eval(parse(text=paste("ifelse(is.list(ml[[",j,"]]),ifelse(!is.null(ml[[",j,"]]$output),ml[[",j,"]]$output,blank),strsplit(ml[[",j,"]],'",xx,"')[[1]][1])",sep=""))))
   }
model.filenames=model.filenames[model.filenames!=""]
purge("inp",model.filenames,ask)
purge("out",model.filenames,ask)
purge("res",model.filenames,ask)
purge("vcv",model.filenames,ask)
purge("rda",model.filenames,ask)
#
#  Delete any tmp files created by export.models or input file creation
#
file.list=list.files(pattern="mark[0123456789]+[:YVX:][:.:]tmp")
unlink(file.list)
file.list=list.files(pattern="markxxx[0123456789][:.:]tmp")
unlink(file.list)
invisible()
}
