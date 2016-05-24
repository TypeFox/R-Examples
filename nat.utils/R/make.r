#' Run a command if input files are newer than outputs
#' 
#' @details \code{cmd} can be an R \code{expression}, which is 
#'   \code{\link{eval}}uated if necessary in the environment calling 
#'   \code{RunCmdForNewerInput}, a string to be passed to \code{\link{system}} 
#'   or NULL/NA in which cases the files are checked and \code{TRUE} or 
#'   \code{FALSE} is returned depending on whether action is required.
#'   
#'   When \code{UseLock=TRUE}, the lock file created is called outfiles[1].lock
#'   
#'   When \code{ReturnInputTimes=TRUE}, the input mtimes are returned as an 
#'   attribute of a logical value (if available).
#'   
#' @param cmd An \code{\link{expression}}, a string or NA/NULL
#' @param infiles Character vector of path to one or more input files
#' @param outfiles Character vector of path to one or more output files
#' @param Verbose Write information to consolse (Default FALSE)
#' @param UseLock Stop other processes working on this task (Default FALSE)
#' @param Force Ignore file modification times and always produce output
#'    if input files exist.
#' @param ReturnInputTimes Return mtimes of input files (default FALSE)
#' @param ... additional parameters passed to \code{\link{system}} call.
#' @return logical indicating if cmd was run or for an R expression, eval(cmd)
#' @export
#' @seealso \code{\link{makelock}, \link{eval}, \link{expression}}
#' @examples \dontrun{
#' RunCmdForNewerInput(expression(myfunc("somefile")))
#' }
RunCmdForNewerInput<-function(cmd,infiles,outfiles,Verbose=FALSE,UseLock=FALSE,
    Force=FALSE, ReturnInputTimes=FALSE,...){
  # note that cmd can be an R expression as in 
  # RunCmdForNewerInput(expression(myfunc("somefile")))
  fisi=NULL
  if(length(infiles)==0){
    if(Verbose) cat("no input files\n")
    return (FALSE)
  } else if(!all(fei<-file.exists(infiles))){
    if(Verbose) cat("some input files missing: ",infiles[!fei],"\n")
    return (FALSE)
  } else if(Force) {
    # do nothing just fall through to end
    if(Verbose) cat("Force=TRUE\n")
  } else if(!all(feo<-file.exists(outfiles))){
    # do nothing just fall through to end
    if(Verbose) cat("outfiles: ",outfiles[!feo],"missing\n")
  } else if( (mit<-max(fisi<-file.info(infiles)$mtime)) <=
             (mot<-min(fiso<-file.info(outfiles)$mtime)) ){
    # check times
    if(Verbose) cat("Skipping",outfiles,"because input files are older\n")
    return(FALSE) 
  } else {
    if(Verbose){
      cat("Overwriting",outfiles,"because 1 or more input files are newer\n")
      cat("Newest input mtime:",strftime(mit),
        "Oldest output mtime:",strftime(mot),"\n")
    } 
  }
  lockfile=paste(outfiles[1],sep=".","lock")
  # return FALSE to signal output doesn't exist
  if(UseLock){
    if(makelock(lockfile))
      on.exit(unlink(lockfile))
    else {
      if(Verbose) cat("Skipping ",outfiles," because someone else is working on ",
        ifelse(length(outfiles)==1,"it","them"),"\n",sep="")
      return(FALSE)
    }
  }
  if(is.expression(cmd)){
    return(eval(cmd, envir=parent.frame()))
  } else if(is.character(cmd)){
    system(cmd,...)
  }
  if(ReturnInputTimes && !is.null(fisi))
    return(structure(TRUE,fisi=fisi))
  else return(TRUE)
}
