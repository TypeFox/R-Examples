RTOOLS <- "Rtools28.exe"
RTOOLS.HOME <- "http://www.murdoch-sutherland.com/Rtools/"
#RTOOLS.URL <- "http://www.murdoch-sutherland.com/Rtools/Rtools28.exe"
RTOOLS.URL <- paste( RTOOLS.HOME, RTOOLS, sep="" )
## Installs to 'C:\Rtools\bin\sh.exe'
SH.EXE <- "C:/Rtools/bin/sh.exe"

permission <- function( msg ) {
  ##library( tcltk )
  res <- tkmessageBox(message=msg,icon="question",type="yesno",default="yes")
  return( tclvalue(res) == "yes" )
}

#isWindows <- getFromNamespace( "isWindows", "pbatR" )
isWindows <- function()
  return( tolower(Sys.info()["sysname"]) == "windows" )
isMac <- function()
  return( tolower(Sys.info()["sysname"]) == "darwin" )
isLinux <- function()
  return( tolower(Sys.info()["sysname"]) == "linux" )

## Returns default install paths
fbat.install.default.exePath <- function() {
  if( isWindows() ) return( "C:/" )
  if( isMac() || isLinux() ) return( "~/bin" )

  stop( "Could not determine system OS type." )
}

## Returns download paths
fbat.install.http <- function() {
  #if( isWindows() ) return( "http://www.people.fas.harvard.edu/~tjhoffm/fbati/fbatWin202c.zip" )
  #if( isLinux() ) return( "http://www.people.fas.harvard.edu/~tjhoffm/fbati/fbatLinux202c.zip" )
  #if( isMac() ) return( "http://www.people.fas.harvard.edu/~tjhoffm/fbati/fbatDarwin202c.zip" )
  if( isWindows() ) return( "http://www.biostat.harvard.edu/%7Efbat/software/fbat202C_win.zip" )
  if( isLinux() ) return( "http://www.biostat.harvard.edu/%7Efbat/software/fbat202C_linux.tar.gz" )
  if( isMac() ) return( "http://www.biostat.harvard.edu/%7Efbat/software/fbat202C_darwin.tar.gz" )

  stop( "Could not determine system OS type." )
}
fbat.install.http.dest <- function() {
  #return("fbat.zip")
  if( isWindows() ) return( "fbat.zip" )
  return( "fbat.tgz" )
}

fbat.exename <- function() {
  #if( isLinux() ) return( "fbatLinux202c" )
  #if( isWindows() ) return( "fbatWin202c.exe" )
  #if( isMac() ) return( "fbatDarwin202c" )

  if( isLinux() || isMac() ) return( "fbat" )
  if( isWindows() ) return( "fbat.exe" )

  stop( "Could not determine system OS type." )
}

fbat.full.exename <- function( exePath=NULL ) {
  fbatname <- fbat.exename()

  if( is.null(exePath) ) exePath <- fbat.install.default.exePath()
  return( paste( exePath, "/", fbatname, sep="" ) )
}

fbat.extrapath <- function() {
  if( isWindows() ) return( "fbat202C_win/" )
  if( isLinux() ) return( "fbat202C_linux/" )
  if( isMac() ) return( "fbat202C/" )

  stop( "Could not determine system OS type." )
}

fbat.install <- function( exePath=NULL ) {
  ## Destination files...
  dest <- fbat.full.exename( exePath )
  if( !file.exists(dest) ) {
    ## Get permission from user.
    if( !permission(paste( "May I download the file [", fbat.install.http(), "], and extract the FBAT executable to the folder [", fbat.full.exename(exePath), "]? You must agree to the FBAT license.", sep="" )) )
      stop( paste( "User did not wish to download the file [", fbat.install.http(), "] to the current working directory. You will need to download and extract fbat manually to the location [", fbat.full.exename(exePath), "].", sep="" ) )

    ## Prep the install directory
    if( is.null( exePath ) ) exePath <- fbat.install.default.exePath()
    if( !file.exists(exePath) ) dir.create( exePath )

    ## Download the file
    cat( "Downloading FBAT...\n" )
    download.file( url=fbat.install.http(), destfile=fbat.install.http.dest(), mode="wb" )
    if( !file.exists( fbat.install.http.dest() ) )
      stop( paste( "Could not download the file [", fbat.install.http(), "] to the current working directory. You may need to download and extract fbat manually to the location [", fbat.full.exename(exePath), "].", sep="" ) )

    ## Extract the file
    cat( "Extracting FBAT...\n" )
    if( !isWindows() ) {
      ## Make it pass codetools check
      assign( "zip.unpack", function() { stop("zip.unpack only available in windows.") } )
    }

    if( isWindows() ) {
      ##zip.file.extract( paste(fbat.extrapath(),fbat.exename(),sep=""), "fbat.zip" )
      ## CAUSES CODETOOLS ERROR, I DO NOT KNOW HOW TO FIX THIS
      zip.unpack( "fbat.zip", "." )
      file.copy( from=paste(fbat.extrapath(),fbat.exename(),sep=""), to="." )
      if( !file.exists( fbat.exename()) )
        stop( paste( "Could not extract the file [fbat.zip] in the current working directory to the file [", fbat.exename(), "]. You may need to download and extract fbat manually to the location [", fbat.full.exename(exePath), "].", sep="" ) )
    }else{
      cat( "tar -zxvf fbat.tgz", "\n" )
      system( "tar -zxvf fbat.tgz" )
      cat( paste( "mv ", fbat.extrapath(), "* .", sep="" ), "\n" )
      system( paste( "mv ", fbat.extrapath(), "* .", sep="" ) )
      if( !file.exists(fbat.exename()) )
        stop( paste( "Could not extract the file [fbat.tgz] in the current working directory to the file [", fbat.exename(), "]. You may need to download and extract fbat manually to the location [", fbat.full.exename(exePath), "].", sep="" ) )
    }

    ## And move it to the proper location
   file.copy( fbat.exename(), dest )
   if( !isWindows() )
    Sys.chmod( dest ) ## OK, this is just stupid. Thank you.s

    ## Make sure it exists then...
    if( !file.exists( dest ) )
      stop( "Failed to install fbat" )
  }
}


sh.install <- function() {
  if( !isWindows() || file.exists(SH.EXE) )
    return()

  #if( !permission( paste( "Can I download and then run the installer for Rtools? This is required for the sh.exe program to control the FBAT program. I would like to download the file [", RTOOLS.URL, "], run it, you will proceed through the installation, and return here. Do not change the defaults, as then this will no longer work. Proceed?", sep="" ) )
  if( !permission( paste( "Can I download and then run the installer for Rtools? This is required for the sh.exe program to control FBAT. This would involve downloading [", RTOOLS.URL, "], running it, installing to the default location, and then returning here. Note that a newer version may be available on the website if this fails, or you prefer that version. Proceed?", sep="" ) ) )
    stop( "Could not install the 'sh.exe' file, terminating." )

  cat( "Downloading file...\n" )
  download.file( url=RTOOLS.URL, destfile="rtools.exe", mode="wb" )
  cat( "Installing file...\n" )
  system( "rtools.exe" )

  if( !file.exists(SH.EXE) )
    stop( "Installation of sh.exe failed. Installation failed." )
}


## DEBUG ONLY
fbatcinstall.debug <- function() {
  ## Do we recognize the OS?
  cat( "isWindows", isWindows(), "\n" )
  cat( "isLinux", isLinux(), "\n" )
  cat( "isMac", isMac(), "\n" )

  ## What's the exename?
  cat( "exename:", fbat.exename(), "\n" )
  cat( "full exename (default):", fbat.full.exename(), "\n" )

  ## Try installing fbat
  fbat.install()
  print( dir() ) ## List directory contents to see what happened.

  ## Try installing sh (windows only)
  sh.install()
}
#fbatcinstall.debug()
