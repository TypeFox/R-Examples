
#.onLoad hook for when namespace is loaded
.onLoad <- function(libname, pkgname) {
  configureExifTool(libname, pkgname, quiet=FALSE, forceInstall=FALSE)
}

configureExifTool <- function(libname, pkgname, quiet=TRUE, forceInstall=FALSE) {
  #find command, will throw error if perl is not installed
  command <- findExifToolCommand(libname, pkgname, quiet=quiet, forceInstall=forceInstall)
  if(is.null(command) || forceInstall) {
    #try downloading/installing exiftool
    if(forceInstall && !quiet) message("Forcing download/install")
    #always message before download/install
    message("ExifTool not found, attempting to install from paleolimbot.github.io/exifr")
    installExifTool(libname, pkgname, quiet=quiet)
    command <- findExifToolCommand(libname, pkgname, quiet=quiet, forceInstall=forceInstall)
  }

  if(is.null(command)) {
    stop("Could not find ExifTool: installation may have failed")
  } else {
    options(exifr.exiftoolcommand=command)
  }
}

findExifToolCommand <- function(libname, pkgname, quiet=TRUE, forceInstall=FALSE) {
  #try straight-up exiftool command
  if(!quiet) message("Trying 'exiftool' on the console...")
  if(testCommand("exiftool --version") && !forceInstall) {
    if(!quiet) message("Found")
    return("exiftool")
  } else {
    #find Perl, check default install location
    if(!quiet) message("Checking Perl installation")
    perlpaths <- c(options("exifr.perlpath")$exifr.perlpath, "perl",
                   "C:\\Perl64\\bin\\perl", "C:\\Perl\\bin\\perl",
                   "c:\\Strawberry\\perl\\bin\\perl")
    perlpath <- NULL
    for(pp in perlpaths) {
      if(testCommand(paste(pp, "--version"))) {
        #found perl
        perlpath <- pp
        break
      }
    }
    if(is.null(perlpath)) {
      stop("Perl is required for exifr and was not found at any of the following locations: ",
           paste(perlpaths, collapse=" "),
           ". Specify perl location by setting options(exifr.perlpath='my/path/to/perl')")
    }
    if(!quiet) message("Found perl at ", perlpath, "; looking for ExifTool")
    #look for exiftool/exiftool.pl in two locations:
    #file.path(libname, pkgname, "exiftool/exiftool.pl"), and file.path("~", "exiftool/exiftool.pl")
    #remember to shellquote filenames!
    commands <- c(paste(perlpath, shQuote(file.path(libname, pkgname, "exiftool/exiftool.pl"))),
                  paste(perlpath, shQuote(file.path("~", "exiftool/exiftool.pl"))))

    for(command in commands) {
      if(testCommand(paste(command, "--version"))) {
        if(!quiet) message("Found valid exiftool command: ", command)
        return(command)
      }
    }

    #perl is installed, but not ExifTool
    if(!quiet) message("No valid ExifTool installation")
    return(NULL)
  }
}

testCommand <- function(command) {
  suppressWarnings(suppressMessages(0==try(system(command, ignore.stdout = TRUE,
                                                  ignore.stderr = TRUE, show.output.on.console = FALSE), silent=TRUE)))
}

installExifTool <- function(libname, pkgname, quiet=TRUE) {
  installlocs <- c(file.path(libname, pkgname), path.expand("~"))
  exiftoolurl <- "http://paleolimbot.github.io/exifr/exiftool.zip"
  #check for a writable install location
  if(!quiet) message("Checking for writable install location...")
  installloc <- NULL
  for(il in installlocs) {
    testfile <- file.path(il, ".dummyfile")
    file.create(testfile, showWarnings = FALSE)
    if(file.exists(testfile)) {
      unlink(testfile)
      installloc <- il
      break
    }
  }
  installloc <- installlocs[1]

  if(is.null(installloc)) {
    stop("Could not find writable install location to install ExifTool: tried ",
         paste(installlocs, collapse=" "))
  }

  if(!quiet) message("Attempting to install ExifTool to ", installloc)

  #download exiftool zip
  tryCatch({
    zipfile <- file.path(installloc, "exiftool.zip")
    if(!quiet) message("Downloading ExifTool from ", exiftoolurl)
    utils::download.file(exiftoolurl, zipfile, quiet = TRUE)
    #check download
    if(!file.exists(zipfile)) {
      stop("Error downloading ExifTool from ", exiftoolurl)
    }
    #unzip
    if(!quiet) message("Extracting ExifTool to ", installloc)
    utils::unzip(zipfile, exdir=installloc, overwrite = TRUE)
    #remove zip file
    unlink(zipfile)
    #return success
    return(TRUE)
  }, error=function(e){
    stop("Error installing ExifTool: ", e)
  })

}
