# save and save.image for ff
# (c) 2009 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2009-10-25
# Last changed: 2009-10-25

# source("d:/mwp/eanalysis/ff/R/ffsave.R")

if (FALSE){
  library(ff)
  file.remove("d:/tmp/a.ff")
  file.remove("d:/tmp/b.ff")
  file.remove("d:/tmp/x.ff")
  file.remove("d:/tmp/y.ff")
  file.remove("d:/tmp/z.ff")

  file.remove("d:/tmp/x.RData")
  file.remove("d:/tmp/x.ffData")
  file.remove("d:/tmp/y.RData")
  file.remove("d:/tmp/y.ffData")
  file.remove("d:/tmp/z.RData")
  file.remove("d:/tmp/z.ffData")

  message("let's create some ff objects")
  n <- 8e3
  a <- ff(sample(n, n, TRUE), vmode="integer", length=n, filename="d:/tmp/a.ff")
  b <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/b.ff")
  x <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/x.ff")
  y <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/y.ff")
  z <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/z.ff")
  df <- ffdf(x=x, y=y, z=z)
  rm(x,y,z)

  message("save all of them")
  ffsave.image("d:/tmp/x")
  str(ffinfo("d:/tmp/x"))

  message("save some of them with shorter relative pathnames ...")
  ffsave(a, b, file="d:/tmp/y", rootpath="d:/tmp")
  str(ffinfo("d:/tmp/y"))

  message("... and add others later")
  ffsave(df, add=TRUE, file="d:/tmp/y", rootpath="d:/tmp")
  str(ffinfo("d:/tmp/y"))

  message("... and add others later")
  system.time(ffsave(a, file="d:/tmp/z", move=TRUE))
  ffinfo("d:/tmp/z")

  message("let's delete/close/remove all objects")
  close(a)  # no file anymore, since we moved a into the ffarchive
  delete(b, df)
  rm(df, a, b, n)
  message("prove it")
  ls()

  message("restore all but ff files in a different directory")
  system.time(ffload("d:/tmp/x", rootpath="d:/tmp2"))
  lapply(ls(), function(i)filename(get(i)))

  delete(a, b, df)
  rm(df, a, b)

  ffdrop(c("d:/tmp/x", "d:/tmp/y", "d:/tmp/z"))
}



#! \name{ffsave}
#! \alias{ffsave}
#! \alias{ffsave.image}
#! \title{
#! Save R and ff objects
#! }
#! \description{
#!  \code{ffsave} writes an external representation of R and ff objects to an \code{ffarchive}.
#!  The objects can be read back from the file at a later date by using the function \code{\link{ffload}}.
#! }
#! \usage{
#! ffsave(...
#! , list = character(0L)
#! , file = stop("'file' must be specified")
#! , envir = parent.frame()
#! , rootpath = NULL
#! , add = FALSE
#! %, overwrite = FALSE
#! , move = FALSE
#! , compress = !move
#! , compression_level = 6
#! , precheck=TRUE
#! )
#! ffsave.image(file = stop("'file' must be specified"), safe = TRUE, ...)
#! }
#! \arguments{
#!   \item{\dots}{
#!   For \code{ffsave} the names of the objects to be saved (as symbols or character strings),
#!   for \code{ffsave.image} further arguments passed to \code{ffsave}
#! }
#!   \item{list}{
#!   A character vector containing the names of objects to be saved.
#! }
#!   \item{file}{
#!   A name for the the \code{ffarchive}, i.e. the two files \code{<file>.RData} and \code{<file>.ffData}
#! }
#!   \item{envir}{
#!   environment to search for objects to be saved.
#! }
#!   \item{add}{
#!   logical indicating whether the objects shall be added to the \code{ffarchive} (in this case \code{rootpath} is taken from an existing archive)
#! }
#! %  \item{overwrite}{
#! %  logical indicating whether an existing archive may be overwritten
#! %}
#!   \item{move}{
#!   logical indicating whether ff files shall be moved instead of copied into the \code{<file>.ffData}
#! }
#!   \item{compress}{
#!   logical specifying whether saving to a named file is to use compression.
#! }
#!   \item{compression_level}{
#!   compression level passed to \code{zip}, default 6
#! }
#!   \item{rootpath}{
#!   optional path component that all \emph{all} ff files share and that can be dropped/replaced when calling \code{\link{ffload}}
#! }
#!   \item{precheck}{
#!  logical: should the existence of the objects be checked before starting to save (and in particular before opening the file/connection)?
#! }
#!   \item{safe}{
#!  logical. If \code{TRUE}, a temporary file is used for creating the saved workspace. The temporary file is renamed to \code{<file>.ffData} if the save succeeds.
#!  This preserves an existing workspace \code{<file>.ffData} if the save fails,
#!  but at the cost of using extra disk space during the save.
#! }
#! }
#! \details{
#!   \code{ffsave} stores objects and ff files in an \code{ffarchive} named \code{<file>}:
#!   i.e. it saves all specified objects via \code{\link{save}} in a file named \code{<file>.RData}
#!   and saves all ff files related to these objects in a zipfile named \code{<file>.ffData} using an external \code{zip} utility.
#!   \cr
#!   By default files are stored relative to the \code{rootpath="\"} and will be restored relative to \code{"\"} (in its original location).
#!   By providing a partial path prefix via argument \code{rootpath} the files are stored relative to this \code{rootpath}.
#!   The \code{rootpath} is stored in the \code{<file>.RData} with the name \code{.ff.rootpath}.
#!   I.e. even if the ff objects were saved with argument \code{rootpath} to \code{ffsave},
#!   \code{\link{ffload}} by default restores in the original location.
#!   By using argument \code{rootpath} to \code{ffload} you can restore relative to a different \code{rootpath}
#!   (and using argument \code{rootpath} to \code{ffsave} gave you shorter relative paths)
#!   \cr
#!   By using argument \code{add} in \code{ffsave} you can add more objects to an existing \code{ffarchive}
#!   and by using argument \code{list} in \code{ffload} you can selectively restore objects.
#!   \cr
#!   The content of the \code{ffarchive} can be inspected using \code{\link{ffinfo}} before actually loading any of the objects.
#!   \cr
#!   The \code{ffarchive} can be deleted from disk using \code{\link{ffdrop}}.
#! }
#! \value{
#!   a character vector with messages returned from the \code{zip} utility (one for each ff file zipped)
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \note{
#!   The ff files are not platform-independent with regard to byte order.
#!   For large files and the zip64 format use \code{zip 3.0} and \code{unzip 6.0} from \url{http://www.info-zip.org/}.
#! }
#! \seealso{
#!   \code{\link{ffinfo}} for inspecting the content of the \code{ffarchive} \cr
#!   \code{\link{ffload}} for loading all or some of the \code{ffarchive} \cr
#!   \code{\link{ffdrop}} for deleting one or more \code{ffarchives}
#! }
#! \examples{
#!   \dontrun{
#!   message("let's create some ff objects")
#!   n <- 8e3
#!   a <- ff(sample(n, n, TRUE), vmode="integer", length=n, filename="d:/tmp/a.ff")
#!   b <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/b.ff")
#!   x <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/x.ff")
#!   y <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/y.ff")
#!   z <- ff(sample(255, n, TRUE), vmode="ubyte", length=n, filename="d:/tmp/z.ff")
#!   df <- ffdf(x=x, y=y, z=z)
#!   rm(x,y,z)
#!
#!   message("save all of them")
#!   ffsave.image("d:/tmp/x")
#!   str(ffinfo("d:/tmp/x"))
#!
#!   message("save some of them with shorter relative pathnames ...")
#!   ffsave(a, b, file="d:/tmp/y", rootpath="d:/tmp")
#!   str(ffinfo("d:/tmp/y"))
#!
#!   message("... and add others later")
#!   ffsave(df, add=TRUE, file="d:/tmp/y", rootpath="d:/tmp")
#!   str(ffinfo("d:/tmp/y"))
#!
#!   message("... and add others later")
#!   system.time(ffsave(a, file="d:/tmp/z", move=TRUE))
#!   ffinfo("d:/tmp/z")
#!
#!   message("let's delete/close/remove all objects")
#!   close(a)  # no file anymore, since we moved a into the ffarchive
#!   delete(b, df)
#!   rm(df, a, b, n)
#!   message("prove it")
#!   ls()
#!
#!   message("restore all but ff files in a different directory")
#!   system.time(ffload("d:/tmp/x", rootpath="d:/tmp2"))
#!   lapply(ls(), function(i)filename(get(i)))
#!
#!   delete(a, b, df)
#!   rm(df, a, b)
#!
#!   ffdrop(c("d:/tmp/x", "d:/tmp/y", "d:/tmp/z"))
#!   }
#! }
#! \keyword{ IO }
#! \keyword{file}

ffsave <-
function (
  ...
, list = character(0L)
, file = stop("'file' must be specified")
, envir = parent.frame()
, rootpath = NULL  # uses getwd()
, add = FALSE
, move = FALSE
, compress = !move
, compression_level = 6
, precheck = TRUE
#, overwrite = FALSE
)
{
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")

    opts <- getOption("save.defaults")
    if (missing(compress) && !is.null(opts$compress))
        compress <- opts$compress
    if (missing(compression_level) && !is.null(opts$compression_level))
        compression_level <- opts$compression_level

    cwd <- getwd()
    on.exit(setwd(cwd))

    # make sure the *.RData and *.ffData files have absolute paths
    dfile <- dirname(file)
    bfile <- basename(file)
    if (!file.exists(dfile))
      dir.create(dfile, recursive=TRUE)
    setwd(dfile)
    dfile <- getwd()
    setwd(cwd)  # looks silly but prevents problems with upper/lower case
    file <- file.path(dfile, bfile)
    # file problem with file.path
    file <- gsub("/+","/",file)
    zipfile <- paste(file, "ffData", sep=".")
    imgfile <- paste(file, "RData", sep=".")

    names <- as.character(substitute(list(...)))[-1L]
    list <- c(list, names)
    if (precheck) {
        ok <- unlist(lapply(list, exists, envir = envir))
        if (!all(ok)) {
            n <- sum(!ok)
            stop(sprintf(ngettext(n, "object %s not found",
              "objects %s not found"), paste(sQuote(list[!ok]),
              collapse = ", ")), domain = NA)
        }
    }

    # filter ffdf and ff only
    filelist <- unlist(lapply(list, function(i){
      x <- get(i, envir = envir)
      if (is.ffdf(x)){
        io <- is.open(x)
        if (is.na(io) || io)
          close(x)
        unlist(lapply(physical(x), filename))
      }else if(is.ff(x)){
        if (is.open(x))
          close(x)
        filename(x)
      }
    }))
    # filter out duplicates
    filelist <- unique(filelist)
    
    if (add){
      tempenvir <- new.env()
      oldlist<- load(imgfile, tempenvir)
      rootpath <- get(".ff.rootpath", tempenvir)
      setwd(rootpath)
      save(list=c(oldlist, list), envir=tempenvir, file=imgfile, compress=compress, precheck=FALSE)
    }else{
      # make sure we have a rootpath and it has absolute path
      if (is.null(rootpath))
        rootpath <- "/"
      setwd(dirname(filelist[[1]]))
      setwd(rootpath)
      rootpath <- sub("//$","/",paste(getwd(), "/", sep=""))

      #if (!overwrite && file.exists(imgfile))
      #  stop("must not overwrite '", imgfile, "'")
      if (file.exists(zipfile)){
        #if (overwrite)
        #  stop("must not overwrite '", zipfile, "'")
        #else
          file.remove(zipfile)
      }

      assign(".ff.rootpath", rootpath, envir=envir)
      on.exit(rm(list=".ff.rootpath", envir=envir), add=TRUE)
      savecall <- match.call(save)
      savecall[[1]] <- as.name("save")
      savecall$file <- imgfile
      savecall$list <- c(".ff.rootpath", list)
      savecall$move <- NULL
      savecall$rootpath <- NULL
      savecall$precheck <- FALSE
      savecall$add <- NULL
      savecall$move <- NULL
      eval(savecall, envir=envir)
    }

    # check rootpath compatibility
    i <- grepl(paste("^", rootpath, sep=""), filelist)
    if (!all(i)){
      print(filelist[!i])
      stop("the previous files do not match the rootpath (case sensitive)")
    }
    nlist <- nchar(filelist)
    nroot <- nchar(rootpath)
    filelist <- substr(filelist, nroot+1, nlist)

    cmd <- paste('zip -@ -', if (compress) compression_level  else 0, if (move) ' -m', ' "', zipfile, '"', sep="")
    ret <- system(cmd, input=filelist, intern=TRUE)
    ret
}


ffsave.image <-
function (
  file = stop("'file' must be specified")
, safe = TRUE
, ...
)
{
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")

    opts <- getOption("save.image.defaults")
    if (is.null(opts))
        opts <- getOption("save.defaults")
    if (missing(safe) && !is.null(opts$safe))
        safe <- opts$safe

    cwd <- getwd()
    on.exit(setwd(cwd))
    # make sure the *.RData and *.ffData files have absolute paths
    dfile <- dirname(file)
    bfile <- basename(file)
    if (!file.exists(dfile))
      dir.create(dfile, recursive=TRUE)
    setwd(dfile)
    dfile <- getwd()
    file <- file.path(dfile, bfile)
    # file problem with file.path
    file <- gsub("/+","/",file)
    setwd(cwd)  # looks silly but prevents problems with upper/lower case
    on.exit()

    if (safe) {
        outfile <- paste(file, "Tmp", sep = "")
        imgfile <- paste(outfile, "RData", sep=".")
        zipfile <- paste(outfile, "ffData", sep=".")
        i <- 0
        while (file.exists(imgfile) || file.exists(zipfile) ) {
            i <- i + 1
            outfile <- paste(file, "Tmp", i, sep = "")
            imgfile <- paste(outfile, "RData", sep=".")
            zipfile <- paste(outfile, "ffData", sep=".")
        }
    }
    else{
      outfile <- file
      imgfile <- paste(outfile, "RData", sep=".")
      zipfile <- paste(outfile, "ffData", sep=".")
    }
    on.exit(file.remove(c(imgfile, zipfile)))

    ret <- ffsave(list = ls(envir = .GlobalEnv, all.names = TRUE)
    , file = outfile
    , envir = .GlobalEnv
    , ...
    )

    on.exit()
    if (safe){
        if (!file.move(imgfile, paste(file, "RData", sep="."))) {
            stop("R image could not be renamed and is left in ", imgfile)
        }
        if (!file.move(zipfile, paste(file, "ffData", sep="."))) {
            stop("ff image could not be renamed and is left in ", zipfile)
        }
    }
    ret
}


#! \name{ffinfo}
#! \alias{ffinfo}
#! \title{
#!   Inspect content of ff saves
#! }
#! \description{
#!   Find out which objects and ff files are in a pair of files saved with \code{\link{ffsave}}
#! }
#! \usage{
#! ffinfo(file)
#! }
#! \arguments{
#!   \item{file}{
#!   a character string giving the name (without extension) of the \code{.RData} and \code{.ffData} files to load
#! }
#! }
#! \value{
#!   a list with components
#!   \item{RData}{a list, one element for each object (named like the object): a character vector with the names of the ff files }
#!   \item{ffData}{a list, one element for each path (names like the path): a character vector with the names of the ff files }
#!   \item{rootpath}{ the root path relative to which the files are stored in the .ffData zip }
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \note{
#!   For large files and the zip64 format use \code{zip 3.0} and \code{unzip 6.0} from \url{http://www.info-zip.org/}.
#! }
#! \seealso{
#!   \code{\link{ffsave}}, \code{\link{ffload}}, \code{\link{ffdrop}}
#! }
#! \keyword{ IO }
#! \keyword{file}

ffinfo <-
function (file)
{
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")

    imgfile <- paste(file, "RData", sep=".")
    zipfile <- paste(file, "ffData", sep=".")
    #if (.Platform$OS.type == "windows"){
    #  zipfile <- chartr("/", "\\", zipfile)
    #}

    envir <- new.env()
    load(imgfile, envir=envir)

    nam <- ls(all.names=TRUE, envir=envir)
    nam <- nam[nam!=".ff.rootpath"]

    rootpath <- get(".ff.rootpath", envir=envir)

    RData <- vector("list", length(nam))
    names(RData) <- nam

    for (i in nam){
      x <- get(i, envir = envir)
      if (is.ffdf(x)){
        RData[[i]] <- unlist(lapply(physical(x), filename))
      }else if(is.ff(x)){
        RData[[i]] <- filename(x)
      }else
        RData[[i]] <- NULL
    }

    cmd <- paste('unzip -Z -1 "', zipfile, '"', sep="")
    ffData <- system(cmd, intern=TRUE)
    spf <- splitPathFile(ffData)
    list(RData=RData, ffData=split(spf$file, spf$path), rootpath=rootpath)
}



#! \name{ffload}
#! \alias{ffload}
#! \title{
#!  Reload ffSaved Datasets
#! }
#! \description{
#!   Reload datasets written with the function \code{ffsave} or \code{ffsave.image} .
#! }
#! \usage{
#! ffload(file, list = character(0L), envir = parent.frame()
#! , rootpath = NULL, overwrite = FALSE)
#! }
#! \arguments{
#!   \item{file}{
#!   a character string giving the name (without extension) of the \code{.RData} and \code{.ffData} files to load
#! }
#!   \item{list}{
#!   An optional vector of names selecting those objects to be restored (default NULL restores all)
#! }
#!   \item{envir}{
#!   the environment where the data should be loaded.
#! }
#!   \item{rootpath}{
#!   an optional rootpath where to restore the ff files (default NULL restores in the original location)
#! }
#!   \item{overwrite}{
#!   logical indicating whether possibly existing ff files shall be overwritten
#! }
#! }
#! \details{
#!   \code{\link{ffinfo}} can be used to inspect the contents an ffsaved pair of \code{.RData} and \code{.ffData} files.
#!   Argument \code{list} can then be used to restore only part of the ffsave.
#! }
#! \value{
#!   A character vector with the names of the restored ff files
#! }
#! \note{
#!   The ff files are not platform-independent with regard to byte order.
#!   For large files and the zip64 format use \code{zip 3.0} and \code{unzip 6.0} from \url{http://www.info-zip.org/}.
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \seealso{
#!   \code{\link{load}}, \code{\link{ffsave}}, \code{\link{ffinfo}}, \code{\link{ffdrop}}
#! }
#! \keyword{ IO }
#! \keyword{file}

ffload <-
function (
  file
, list = character(0L)
, envir = parent.frame()
, rootpath = NULL
, overwrite = FALSE
)
{
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")

    cwd <- getwd()
    on.exit(setwd(cwd))

    # make sure the *.RData and *.ffData files have absolute paths
    dfile <- dirname(file)
    bfile <- basename(file)
    setwd(dfile)
    dfile <- getwd()
    setwd(cwd)  # looks silly but prevents problems with upper/lower case
    file <- file.path(dfile, bfile)
    # file problem with file.path
    file <- gsub("/+","/",file)
    zipfile <- paste(file, "ffData", sep=".")
    imgfile <- paste(file, "RData", sep=".")

    cmd <- paste('unzip -Z -1 "', zipfile, '"', sep="")
    ffData <- system(cmd, intern=TRUE)

    haslist <- !!length(list)
    if (!overwrite || haslist || !is.null(rootpath)){
      tempenvir <- new.env()
      load(imgfile, envir=tempenvir)
      oldrootpath <- get(".ff.rootpath", envir=tempenvir)
      if (length(list)){
        isinenv <- sapply(list, exists, envir=tempenvir)
        if (!all(isinenv))
          stop("not in ffarchive: ", paste('"', list[!isinenv]), '"', collapse=",", sep="")
      }else{
        list <- ls(all.names=TRUE, envir=tempenvir)
        list <- list[list!=".ff.rootpath"]
      }
      names(list) <- list
      if (is.null(rootpath)){
        rootpath <- oldrootpath
        list <- unlist(lapply(list, function(i){
          if (!overwrite && exists(i, envir)){
            warning("did not overwrite object '", i, "'")
            ret <- NULL
          }else{
            x <- get(i, envir = tempenvir)
            assign(i, x, envir=envir)
            if (is.ffdf(x)){
              ret <- unlist(lapply(physical(x), filename))
            }else if(is.ff(x)){
              ret <- filename(x)
            }else{
              ret <- NULL
            }
          }
          rm(list=i, envir = tempenvir)
          ret
        }))
      }else{
        if (!file.exists(rootpath))
          dir.create(rootpath, recursive=TRUE)
        #if (!file.exists(oldrootpath))
        #  dir.create(oldrootpath, recursive=TRUE)
        #setwd(oldrootpath)
        #oldrootpath <- getwd()
        #setwd(cwd) # looks silly but prevents problems with upper/lower case
        setwd(rootpath)
        rootpath <- getwd()

        oldrootpathsep <- paste(sub("/$", "", oldrootpath), "/", sep="")
        rootpathsep <- paste(sub("/$", "", rootpath), "/", sep="")

        list <- unlist(lapply(list, function(i){
          if (!overwrite && exists(i, envir)){
            warning("did not overwrite object '", i, "'")
            ret <- NULL
          }else{
            x <- get(i, envir = tempenvir)
            if (is.ffdf(x)){
              ret <- unlist(lapply(physical(x), function(y){
                newnam <- sub(oldrootpathsep, rootpathsep, filename(y))
                physical(y)$filename <- newnam
                newnam
              }))
              assign(i, x, envir=envir)
              ret
            }else if(is.ff(x)){
              newnam <- sub(oldrootpathsep, rootpathsep, filename(x))
              physical(x)$filename <- newnam
              assign(i, x, envir=envir)
              ret <- newnam
            }else{
              assign(i, x, envir=envir)
              ret <- NULL
            }
          }
          rm(list=i, envir = tempenvir)
          ret
        }))
      }
      list <- unique(list)

      if (!length(list))
        return(character())

      # determine those filenames in zip file (e.g. /tmp/x) that match the filename in the ff object (e.g. d:/tmp/x)
      i <- unlist(lapply(ffData, function(x){
        i <- grep(x, list)
        n <- length(i)
        if (n==1){
          if (n>1)
            stop("zip ff name '", x, "'matches multiply in list")
          i
        }else{
          NA
        }
      }))
      j <- match(1:length(list), i)
      if (any(is.na(j)))
        stop("could not match list in zip: ", paste('"', list[is.na(j)]), '"', collapse=",", sep="")

      ffData <- ffData[!is.na(i)]

      if (overwrite)
        cmd <- paste('unzip -o "', zipfile, '" ', paste('"',ffData,'"', sep="", collapse=" "), ' -d "',  rootpath, '"', sep="")
      else
        cmd <- paste('unzip -n "', zipfile, '" ', paste('"',ffData,'"', sep="", collapse=" "), ' -d "', rootpath, '"', sep="")
      ext <- system(cmd, intern=TRUE)[-1]

    }else{

      load(imgfile, envir=envir)

      rootpath <- get(".ff.rootpath", envir=envir)
      rm(".ff.rootpath", envir=envir)

      if (overwrite)
        cmd <- paste('unzip -o -d "',  rootpath ,'" "', zipfile, '"', sep="")
      else
        cmd <- paste('unzip -n -d "', rootpath ,'" "', zipfile, '"', sep="")
      ext <- system(cmd, intern=TRUE)[-1]
    }

    i <- unlist(lapply(ffData, function(x){
      i <- grep(x, ext)
      n <- length(i)
      if (n){
        if (n>1)
          warning("multiple matches of '", x, "' in unzip output")
        TRUE
      }else{
        if (overwrite)
          warning("ERROR: did not extract file '", x, "'")
        else
          warning("NOTE: did not overwrite file '", x, "'")
        FALSE
      }
    }))
    ffData[i]
}



#! \name{ffdrop}
#! \alias{ffdrop}
#! \title{
#!   Delete an ffarchive
#! }
#! \description{
#!   Delete the \code{<file>.Rdata} and \code{<file>.ffData} files behind an \code{ffarchive}
#! }
#! \usage{
#! ffdrop(file)
#! }
#! \arguments{
#!   \item{file}{
#!   vector of archive filenames (without extensions)
#! }
#! }
#! \value{
#!   A list with components
#!   \item{RData}{vector with results of \code{\link{file.remove}} on RData files }
#!   \item{ffData}{Description of 'comp2'}
#! }
#! \author{
#!   Jens Oehlschlägel
#! }
#! \note{
#!   This deletes file on disk without warning
#! }
#! \seealso{
#!   \code{\link{ffsave}}, \code{\link{ffinfo}}, \code{\link{ffload}}
#! }
#! \keyword{ IO }
#! \keyword{file}

ffdrop <- function(file){
    if (!is.character(file) || file == "")
        stop("'file' must be non-empty string")
    zipfile <- paste(file, "ffData", sep=".")
    imgfile <- paste(file, "RData", sep=".")
    RData <- file.remove(imgfile)
    ffData <- file.remove(zipfile)
    names(RData) <- imgfile
    names(ffData) <- zipfile
    list(
       RData=RData
    , ffData=ffData
    )
}
