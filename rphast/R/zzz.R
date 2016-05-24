#.onAttach <- function(libname, pkgname) {
#  if(Sys.info()[[1]] == "Windows") { 
#    regDir <- paste(shell("echo %appdata%", intern=TRUE), "\\rphast\\", sep="")
#  }  else {
#    regDir <- paste(system("echo $HOME", intern=TRUE), "/.rphast/", sep="")
#  }
#  
#  Sys.setenv(rphastRegDir=regDir)
#  
#  if((!file.exists(paste(regDir, "registered", sep=""))) &&
#     (!file.exists(paste(regDir, "notRegistered", sep="")))) {
#    packageStartupMessage("
#************************************************************************
#** If you find RPHAST useful, we would appreciate if you let us know  **
#**   so we can better understand our user base.  Registration is free **
#**   and can be anonymous.                                            **
#** See ?register.rphast for more details, or type nothanks.rphast()   **
#**   to stop these reminders without registering.  Thanks!            **
#************************************************************************
#")
#  }
#}

##' If you are making use of RPHAST, we would appreciate if you would let us
##' know.  This will send your name, email, and institution (all optional),
##' as well as your IP address and rphast version to our server.  We will
##' not share your information with anyone.
##' If you choose to send your email address, we may use it (very rarely) to
##' let you know about major new releases and bug fixes.
##'
##' Once you register, an empty file called "registered" will be created
##' in ~/.rphast (non-Windows) or \%appData\%\\rphast (Windows) which will
##' indicate to us that you have registered, and you will no longer receive
##' any reminders to register when rphast is loaded.
##' @title Register RPHAST
##' @param name Your Name (Optional).  Let us know who you are, if you want.
##' @param email Your Email Address (Optional).  If given, we may use it
##' very rarely to let you know about major new releases and bug fixes.  We
##' will not share your email address with anyone.
##' @param institution Your Institution (Optional).  Let us know where
##' you are from.
##' @param comments Anything else you'd like to tell us!
##' @seealso \code{\link{nothanks.rphast}} to get rid of registration reminders
##' without actually registering.
##' @export
##' @author Nicholas Peterson
register.rphast <- function(name="", email="", institution="", comments="") {
  bName <- URLencode(name, reserved=TRUE);
  bEmail <- URLencode(email, reserved=TRUE);
  bInstitution <- URLencode(institution, reserved=TRUE);
  bComments <- URLencode(comments, reserved=TRUE);
  invisible(read.csv(file=paste("http://compgen.bscb.cornell.edu/rphast/register.php?name=", bName, "&email=", bEmail, "&institution=", bInstitution, "&version=", packageDescription("rphast")$Version, "&isValid=TRUE&comment=", bComments, sep="")))
  regDir <- Sys.getenv("rphastRegDir")
  dir.create(regDir, showWarnings=FALSE)
  invisible(file.create(paste(regDir, "registered", sep=""), showWarnings=FALSE))
}


##' Once this is called rphast will no longer produce startup messages
##' prodding you to register.
##' @title Stop rphast registration reminders
##' @note This creates an empty file called "notRegistered" in
##' \code{Sys.getenv("rphastRegDir")}.  The rphastRegDir is usually in
##' \%appdata\%\\rphast for Windows systems and ~/.rphast for other systems.
##' @export
##' @author Melissa J. Hubisz
nothanks.rphast <- function() {
  cat("You will no longer hear any registration reminders from us!\n")
  regDir <- Sys.getenv("rphastRegDir")
  dir.create(regDir, showWarnings=FALSE)
  invisible(file.create(paste(regDir, "notRegistered", sep=""), showWarnings=FALSE))
}
