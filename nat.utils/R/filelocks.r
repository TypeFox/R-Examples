#' Make and remove (NFS safe) lock files
#'
#' Creates a lock file on disk containing a message that should identify the
#' current R session. Will return FALSE is someone else has already made a
#' lockfile. In order to avoid race conditions typical on NFS mounted drives
#' makelock appends a unique message to the lock file and then reads the file
#' back in. Only if the unique message is the first line in the file will
#' makelock return TRUE.
#' @param lockfile Path to lockfile
#' @param lockmsg Character vector with message to be written to lockfile
#' @param CreateDirectories Recursively create directories implied by lockfile path
#' @return logical indicating success
#' @author jefferis
#' @export
#' @examples
#' makelock(lock<-tempfile())
#' stopifnot(!makelock(lock))
#' removelock(lock)
makelock<-function(lockfile,lockmsg,CreateDirectories=TRUE){
  lockdir=dirname(lockfile)
  if(!file.exists(lockdir)){
    if(CreateDirectories) dir.create(lockdir,recursive=TRUE)
    else stop("Lock Directory for lockfile ",lockfile," does not exist")
  } 
  if(missing(lockmsg)) lockmsg=paste(system('hostname',intern=TRUE),Sys.getenv("R_SESSION_TMPDIR"))
  if (file.exists(lockfile)) return (FALSE)
  # note the use of paste makes the message writing atomic
  cat(paste(lockmsg,"\n",sep=""),file=lockfile,append=TRUE,sep="")
  firstline=readLines(lockfile,n=1)
  if(firstline!=lockmsg){
    # somebody else got there first
    return(FALSE)
  } else return(TRUE)
}

#' Remove lock file
#' 
#' \code{removelock} displays a warning and returns false if lockfile cannot 
#'  be removed. No error message is given if the file does not exist.
#' @export
#' @rdname makelock
removelock<-function(lockfile){
  if(unlink(lockfile)!=0) {
    warning("Unable to remove ",lockfile)
    return (FALSE)
  }
  return (TRUE)
}
