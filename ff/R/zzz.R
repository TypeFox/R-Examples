# Array utilities for ff
# (c) 2007 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2007-08-24
# Last changed: 2007-10-25

# source("d:/mwp/eanalysis/ff/R/zzz.R")

.onLoad <- function(lib, pkg) {
  ##library.dynam("ff", pkg, lib) use useDynLib(ff) in NAMESPACE instead
  ##packageStartupMessage("Loading package ff", packageDescription("ff", lib, fields="Version"), "")
  # allow fftempdir to be set before package is loaded
  if (is.null(getOption("fftempdir"))){
    # create tempdir and make tempdir name independent of platform (otherwise dirname(tempfile)!=getOption("fftempdir"))
    options(fftempdir=standardPathFile(tempdir()))
  }
  if (is.null(getOption("ffextension")))
    options(ffextension="ff")
  if (is.null(getOption("fffinonexit")))
    options(fffinonexit=TRUE)
  if (is.null(getOption("ffpagesize")))
    options(ffpagesize=getdefaultpagesize())
  if (is.null(getOption("ffcaching")))
    options(ffcaching="mmnoflush")
  if (is.null(getOption("ffdrop")))
    options(ffdrop=TRUE)
  if (is.null(getOption("ffbatchbytes"))){
    # memory.limit is windows specific
    if (.Platform$OS.type=="windows")
    {
      if (getRversion()>="2.6.0")  # memory.limit was silently changed from 2.6.0 to return in MB instead of bytes
        options(ffbatchbytes=memory.limit()*(1024^2/100))
      else
        options(ffbatchbytes=memory.limit()/100)
    } else {
      # some magic constant
      options(ffbatchbytes=16*1024^2)
    }
  }
  if (is.null(getOption("ffmaxbytes"))){
    # memory.limit is windows specific
    if (.Platform$OS.type=="windows")
    {
      if (getRversion()>="2.6.0")  # memory.limit was silently changed from 2.6.0 to return in MB instead of bytes
        options(ffmaxbytes=0.5*memory.limit()*(1024^2))
      else
        options(ffmaxbytes=0.5*memory.limit())
    } else {
      # some magic constant
      options(ffmaxbytes=0.5*1024^3)
    }
  }
  # if we want an explicit list of ff objects, we should store them in an environment with hash=TRUE (much faster than a list)
  #assign(".fftemp", new.env(hash=TRUE), envir=globalenv())

  # assign(".vimplemented"
  # , c(
      # boolean   = .Call("ffmode_implemented", .ffmode["boolean"]  , PACKAGE="ff")
    # , logical   = .Call("ffmode_implemented", .ffmode["logical"]  , PACKAGE="ff")
    # , quad      = .Call("ffmode_implemented", .ffmode["quad"]     , PACKAGE="ff")
    # , nibble    = .Call("ffmode_implemented", .ffmode["nibble"]   , PACKAGE="ff")
    # , byte      = .Call("ffmode_implemented", .ffmode["byte"]     , PACKAGE="ff")
    # , ubyte     = .Call("ffmode_implemented", .ffmode["ubyte"]    , PACKAGE="ff")
    # , short     = .Call("ffmode_implemented", .ffmode["short"]    , PACKAGE="ff")
    # , ushort    = .Call("ffmode_implemented", .ffmode["ushort"]   , PACKAGE="ff")
    # , integer   = .Call("ffmode_implemented", .ffmode["integer"]  , PACKAGE="ff")
    # , single    = .Call("ffmode_implemented", .ffmode["single"]   , PACKAGE="ff")
    # , double    = .Call("ffmode_implemented", .ffmode["double"]   , PACKAGE="ff")
    # , complex   = .Call("ffmode_implemented", .ffmode["complex"]  , PACKAGE="ff")
    # , raw       = .Call("ffmode_implemented", .ffmode["raw"]      , PACKAGE="ff")
    # , character = .Call("ffmode_implemented", .ffmode["character"], PACKAGE="ff")
    # ), envir=globalenv()
  # )

  # list of possible coercions without information loss
  # assign(".vcoerceable"
  # , lapply(
      # list(
        # boolean   = 1:13
      # , logical   = c(2L, 5L, 7L, 9:12)
      # , quad      = 3:13
      # , nibble    = 4:13
      # , byte      = c(5L, 7L, 9:12)
      # , ubyte     = 6:13
      # , short     = c(7L, 9:12)
      # , ushort    = 8:12
      # , integer   = 9:12
      # , single    = 10:12
      # , double    = 11:12
      # , complex   = 12L
      # , raw       = 6:12
      # , character = 14L
      # )
      # , function(i)i[.vimplemented[i]]
    # )
    # , envir=globalenv()
  # )
	
  .vimplemented <<- c(
		boolean   = .Call("ffmode_implemented", .ffmode["boolean"]  , PACKAGE="ff")
	, logical   = .Call("ffmode_implemented", .ffmode["logical"]  , PACKAGE="ff")
	, quad      = .Call("ffmode_implemented", .ffmode["quad"]     , PACKAGE="ff")
	, nibble    = .Call("ffmode_implemented", .ffmode["nibble"]   , PACKAGE="ff")
	, byte      = .Call("ffmode_implemented", .ffmode["byte"]     , PACKAGE="ff")
	, ubyte     = .Call("ffmode_implemented", .ffmode["ubyte"]    , PACKAGE="ff")
	, short     = .Call("ffmode_implemented", .ffmode["short"]    , PACKAGE="ff")
	, ushort    = .Call("ffmode_implemented", .ffmode["ushort"]   , PACKAGE="ff")
	, integer   = .Call("ffmode_implemented", .ffmode["integer"]  , PACKAGE="ff")
	, single    = .Call("ffmode_implemented", .ffmode["single"]   , PACKAGE="ff")
	, double    = .Call("ffmode_implemented", .ffmode["double"]   , PACKAGE="ff")
	, complex   = .Call("ffmode_implemented", .ffmode["complex"]  , PACKAGE="ff")
	, raw       = .Call("ffmode_implemented", .ffmode["raw"]      , PACKAGE="ff")
	, character = .Call("ffmode_implemented", .ffmode["character"], PACKAGE="ff")
	)
	
  # list of possible coercions without information loss
  .vcoerceable <<- lapply(
      list(
        boolean   = 1:13
      , logical   = c(2L, 5L, 7L, 9:12)
      , quad      = 3:13
      , nibble    = 4:13
      , byte      = c(5L, 7L, 9:12)
      , ubyte     = 6:13
      , short     = c(7L, 9:12)
      , ushort    = 8:12
      , integer   = 9:12
      , single    = 10:12
      , double    = 11:12
      , complex   = 12L
      , raw       = 6:12
      , character = 14L
      )
      , function(i)i[.vimplemented[i]]
    )
	
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("Attaching package ff")
  packageStartupMessage('- getOption("fftempdir")=="',getOption("fftempdir"),'"\n',sep='')
  packageStartupMessage('- getOption("ffextension")=="',getOption("ffextension"),'"\n',sep='')
  packageStartupMessage('- getOption("ffdrop")==',getOption("ffdrop"),'\n',sep='')
  packageStartupMessage('- getOption("fffinonexit")==',getOption("fffinonexit"),'\n',sep='')
  packageStartupMessage('- getOption("ffpagesize")==',getOption("ffpagesize"),'\n',sep='')
  packageStartupMessage('- getOption("ffcaching")=="',getOption("ffcaching"),'"  -- consider "ffeachflush" if your system stalls on large writes\n',sep='')
  packageStartupMessage('- getOption("ffbatchbytes")==',getOption("ffbatchbytes"),' -- consider a different value for tuning your system\n',sep='')
  packageStartupMessage('- getOption("ffmaxbytes")==',getOption("ffmaxbytes"),' -- consider a different value for tuning your system\n',sep='')
  # if (getRversion()<="2.10.0"){
    # packageStartupMessage('fixing [.AsIs in base namespace because if the NextMethod("[") returns a different class, [.AsIs was reverting this')
    #assignInNamespace(
    #  "[.AsIs"
    #, function (x, i, ...){
    #    ret <- NextMethod("[")
    #    oldClass(ret) <- c("AsIs", oldClass(ret))
    #    ret
    #  }
    #, "base"
    #)
    # assignInNamespace(
      # "[.AsIs"
    # , function (x, i, ...)I(NextMethod("["))
    # , "base"
    # )
  # }
}

.Last.lib <- function(libpath) {
   packageStartupMessage("Detaching package ff")
  # if (getRversion()<="2.10.0"){
    # packageStartupMessage('restoring [.AsIs')
    # assignInNamespace(
      # "[.AsIs"
    # , function (x, i, ...){
        # ret <- NextMethod("[")
        # oldClass(ret) <- c("AsIs", oldClass(x))
        # ret
      # }
    # , "base"
    # )
  # }
}

.onUnload <- function(libpath){
   packageStartupMessage("Unloading package ff")
   #remove(list=".fftemp", envir=globalenv())
   #gc()
   library.dynam.unload("ff", libpath)
   if (unlink(getOption("fftempdir"), recursive = TRUE))
     packageStartupMessage("Error in unlinking fftempdir")
   else
     options(fftempdir=NULL, ffextension=NULL, fffinonexit=NULL, ffpagesize=NULL, ffcaching=NULL, ffdrop=NULL, ffbatchbytes=NULL, ffmaxbytes=NULL)
}
