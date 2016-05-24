"bugs.run" <-
    function(n.burnin, bugs.directory,
             useWINE=.Platform$OS.type != "windows", WINE=NULL,
             newWINE=TRUE, WINEPATH=NULL)
{

  if(useWINE && !is.R())
    stop("Non-Windows platforms not yet supported in R2WinBUGS for S-PLUS")

  ## Is bugs.directory defined in Windows (where second character is :
  ## i.e. C:\Program...) or Unix style path?
  if(useWINE && (substr(bugs.directory, 2, 2) == ":")) {
    bugs.directory <- win2native(bugs.directory, newWINE=newWINE, WINEPATH=WINEPATH)
  }

  ## Update the lengths of the adaptive phases in the Bugs updaters
  try(bugs.update.settings(n.burnin, bugs.directory))

  ## Return the lengths of the adaptive phases to their original settings
  if(is.R()) {
    .fileCopy <- file.copy
  } else {
    .fileCopy <- splus.file.copy
  }
  on.exit(try(.fileCopy(file.path(bugs.directory, "System/Rsrc/Registry_Rsave.odc"),
                        file.path(bugs.directory, "System/Rsrc/Registry.odc"),
                        overwrite=TRUE)))

  ## Search Win*.exe (WinBUGS executable) within bugs.directory
  dos.location <- file.path(bugs.directory,
                            grep("^Win[[:alnum:]]*[.]exe$",
                                 list.files(bugs.directory), value=TRUE)[1])
  if(!file.exists(dos.location))
    stop(paste("WinBUGS executable does not exist in", bugs.directory))

  ## Call Bugs and have it run with script.txt
  bugsCall <- paste("\"", dos.location, "\" /par \"",
                    native2win(file.path(getwd(), "script.txt"),
                               useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH),
                    "\"", sep="")
  if(useWINE) bugsCall <- paste(WINE, bugsCall)
  temp <- system(bugsCall)
  if(temp == -1)
    stop("Error in bugs.run().\nCheck that WinBUGS is in the specified directory.")

  ## Stop and print an error message if Bugs did not run correctly
  if(is.R()) {
    tmp <- scan("coda1.txt", character(), quiet=TRUE, sep="\n")
  } else {
    tmp <- scan("coda1.txt", character(), sep="\n")
  }
  if(length(grep("BUGS did not run correctly", tmp)) > 0)
    stop("Look at the log file and\ntry again with 'debug=TRUE' to figure out what went wrong within Bugs.")
}
