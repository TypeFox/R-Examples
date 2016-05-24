findUnixBinary <- function(x)
{
  ## --- Environmental variable ---
  tmp <- Sys.getenv(toupper(x))
  if(nchar(tmp) != 0 && file.exists(tmp)) return(tmp)
  ## else

  ## --- Standard place ---
  tmp <- paste("/usr/bin", x, sep="")
  if(file.exists(tmp)) return(tmp)
  ## else ...

  ## --- Which ---
  tmp <- system(paste("which ", x, sep=""), intern=TRUE)
  if(length(tmp) != 0 && file.exists(tmp)) return(tmp)
  ## else ..

  ## --- Locate ---
  tmp <- system(paste("locate ", x, " | grep bin/", x, "$", sep=""), intern=TRUE)
  tmp <- tmp[length(tmp)] ## keep only last hit
  if(length(tmp) > 0 && file.exists(tmp)) return(tmp)

  stop(paste("couldn't find", x, "binary file"))
}

native2win <- function(x, useWINE=.Platform$OS.type != "windows",
                       newWINE=TRUE, WINEPATH=NULL)
{
  if (is.R()){
    ## Translate Unix path to Windows (wine) path
    if(useWINE) {
      if(newWINE) {
        if(is.null(WINEPATH)) WINEPATH <- findUnixBinary(x="winepath")
        x <- system(paste(WINEPATH, "-w", x), intern=TRUE)
        gsub("\\\\", "/", x) ## under wine BUGS cannot use \ or \\
      } else {
        winedriveRTr(x)
      }
    } else {
      x
    }
  } else { #S-PLUS
      gsub("\\\\", "/", x)  
  }
}

win2native <- function(x, useWINE=.Platform$OS.type != "windows",
                       newWINE=TRUE, WINEPATH=NULL)
{
  ## Translate Windows path to native (unix) path
  if(useWINE) {
    if(newWINE) {
      if(is.null(WINEPATH)) WINEPATH <- findUnixBinary(x="winepath")
      system(paste(WINEPATH, " \"", x, "\"", sep=""), intern=TRUE)
    } else {
      winedriveTr(x)
    }
  } else {
    x
  }
}

winedriveMap <- function(config="~/.wine/config")
{
  ## Get drive mapping table from ~/.wine/config
  if(!file.exists(config)) return(NULL);
  con <- readLines(config)
  con <- con[- grep("^;", con)]
  drive <- con[grep("^\\[Drive ", con)]
  drive <- substr(drive, 8, 8)
  drive <- paste(drive, ":", sep="")
  path <- con[grep("Path", con)]
  len <- length(drive)
  path <- path[1:len]
  dir <- sapply(path,
                function(x) {
                  foo <- unlist(strsplit(x, "\""))
                  foo[length(foo)]
                })
  dir <- sub("%HOME%",tools::file_path_as_absolute("~"),dir)
  data.frame(drive=I(drive), path=I(dir), row.names=NULL)
}

winedriveTr <- function(windir, DriveTable=winedriveMap())
{
  ## Translate Windows path to native (Unix) path
  win.dr <- substr(windir, 1, 2)
  ind <- pmatch(toupper(win.dr), DriveTable$drive)
  native.dr <- DriveTable$path[ind]
  sub(win.dr, native.dr, windir)
}

winedriveRTr <- function(unixpath, DriveTable=winedriveMap())
{
  ## Translate Unix path to Windows (wine) path
  blocks <- strsplit(unixpath,"/")[[1]]
  cblocks <- c("/",sapply(1+seq(along=blocks[-1]),
                          function(n) paste(blocks[1:n],collapse="/")))
  path <- match(cblocks,DriveTable$path)
  if(any(!is.na(path))) {
    unixdir <- cblocks[which.min(path)]
    windrive <- paste(DriveTable$drive[min(path,na.rm=TRUE)],"/",sep="")
    winpath <- sub("//","/",sub(unixdir,windrive,unixpath)) ## kludge
  } else {
    stop("can't find equivalent Windows path: file may be inaccessible")
  }
  winpath
}
