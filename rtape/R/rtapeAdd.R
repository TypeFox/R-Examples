#' Add object to the tape.
#'
#' This function serializes and appends a given object on the end of the tape file. If the tape file does not exists, it is created.
#'
#' @param fName Name of the tape file.
#' @param what Object to be stored.
#' @param skipNULLs If true and \code{what} is \code{NULL}, nothing is written to the tape. 
#' @param fileFormat File format; should be left default. See \code{\link{guessFileFormat}} and \code{\link{makeFileFormat}} for details.
#' @param safe If \code{"try"} or \code{"retry"}, rtape will use dirlock to ensure that no other rtape safe appending is in progress. In case of conflict, the function in "try" mode immediately returns \code{FALSE} and does not try again, while in "retry" mode it sleeps \code{retryTime} seconds and tries again till the dirlock is opened.
#' @param retryTime If \code{safe} is \code{"retry"}, this parameter sets the interval between writing attempts. Expressed in seconds.
#' @note Remember to use the same \code{fileFormat} value to all writes to a certain tape (or use default format guesser to guarantee this); if not, the tape will become unreadable. For the same reason, don't try to put custom headers/footers or append other data inside tape stream.
#' This function is thread/process safe only if you use \code{safe} mode. However, in this case it may jam on a broken dirlock (for instance when the locking R process crashed during write); you may fix this problem manually by removing the locking dir. Its name is always \code{.rtape_<tape file name>_lock}. Waiting in retry mode is performed via \code{\link{Sys.sleep}}, so it is not a busy wait.
#' @author Miron B. Kursa \email{M.Kursa@@icm.edu.pl}
#' @examples
#' unlink('tmp.tape')
#' #Record something on the tape
#' rtapeAdd('tmp.tape',iris)
#' rtapeAdd('tmp.tape',sin(1:10))
#' #Read whole tape to the list, so we could examine it
#' rtapeAsList('tmp.tape')->stored
#' print(stored)
#' unlink('tmp.tape')

rtapeAdd<-function(fName,what,skipNULLs=FALSE,fileFormat=guessFileFormat(fName),safe=FALSE,retryTime=0.1){
 stopifnot(length(fName)==1)
 lockName<-sprintf('.rtape_%s_lock',fName)
 if(identical(safe,"try")){
  if(!dir.create(lockName,FALSE)) return(invisible(FALSE)) 
   else on.exit({unlink(lockName,TRUE)})
 }else{
  if(identical(safe,"retry")){
    #Block on dirlock
    while(!dir.create(lockName,FALSE)) Sys.sleep(retryTime)
    on.exit({unlink(lockName,TRUE)})
  }
 }
 if(skipNULLs & is.null(what)) return(invisible(NULL))
 fileFormat(fName,open="ab")->con
 serialize(what,con,ascii=FALSE)
 close(con)
 return(invisible(TRUE))
}

