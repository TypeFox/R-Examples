## Modified Function from EmBedBugs package (embedR), by Kenneth Rice
## adapted from http://www.mrc-bsu.cam.ac.uk/personal/ken/embed.html
## old page was: www.statslab.cam.ac.uk/~krice/embed.html
formatData <- function (datalist){
    if (!is.list(datalist) || is.data.frame(datalist)) 
        stop("Argument to format.data must be a list.")
    n <- length(datalist)
    datalist.string <- as.list(rep(NA, n))
    for (i in 1:n) {
        if (length(datalist[[i]]) == 
            1) 
            datalist.string[[i]] <- paste(names(datalist)[i], 
                "=", as.character(datalist[[i]]), sep = "\n")
        if (is.vector(datalist[[i]]) & length(datalist[[i]]) > 
            1) 
            datalist.string[[i]] <- paste(names(datalist)[i], 
                "=c(", paste(as.character(datalist[[i]]), collapse = ", "), 
                ")", sep = "\n")
        if (is.array(datalist[[i]])) 
            datalist.string[[i]] <- paste(names(datalist)[i], 
                "= structure(.Data= c(", paste(as.character(as.vector(aperm(datalist[[i]]))), 
                  collapse = ", "), "), .Dim=c(", paste(as.character(dim(datalist[[i]])), 
                  collapse = ", "), "))", sep = "\n")
    }
    datalist.tofile <- paste("list(", paste(unlist(datalist.string), 
        collapse = ", "), ")", sep = "\n")
    return(datalist.tofile)
}

format4Bugs <- function (dataList, digits=5){
  if (!is.list(dataList) || is.data.frame(dataList)) 
    stop("Argument to formatdata() must be a list.")
#   tmp <- tempfile("dat")
#   on.exit(unlink(tmp))
  ## make sure there is no more than 14 digits
  if (digits > 14) digits <- 14
  dataListString <- lapply(dataList,
                            function(x) {
                              if (is.data.frame(x)) x <- as.matrix(x)
                              dimnames(x) <- NULL
                              if (is.integer(x)) formatC(x, digits=0)
                              ## make sure "E" instead of "e"
                              else formatC(x, format="E", digits=digits)
                            })

  foo <- formatData(dataListString)
  foo <- gsub('"', '', foo)
  foo
# The following method does not always work!  
#   dput(dataListString, tmp)
#   foo <- readLines(tmp)
#   foo <- sub('^structure\\(', '', foo)
#   foo <- gsub('structure\\(c', 'structure\\(.Data=c', foo)
#   foo <- sub('\\), .Names.*$', '\\)', foo)
#   foo <- gsub('"', '', foo)
#   foo <- gsub('as.integer\\(c\\(', 'c\\(', foo)
#   foo <- gsub('\\)\\)', '\\)', foo)
#   foo
}


## generate t.cen from a Surv object

# t.cen <- function(time, status) {
#   ifelse(status == 1, 0, time)
# }

## get drive mapping table from ~/.wine/config
## with changes from Ben Bolker
driveMap <- function(config = "~/.wine/config") {
  if (!file.exists(config)) return (NULL);
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
  dir <- sub("%HOME%",tools:::file_path_as_absolute("~"), dir)
  data.frame(drive = I(drive), path = I(dir), row.names=NULL)
}


## translate windows dir to native dir
driveTr <- function(windir, DriveTable) {
##  .DriveTable <- driveMap(file.path(Sys.getenv("HOME"), ".wine/config"))
##  win.dr <- unlist(strsplit(windir, ":"))[1]
  win.dr <- substr(windir, 1, 2)
  ind <- pmatch(toupper(win.dr), DriveTable$drive)
  native.dr <- DriveTable$path[ind]
  sub(win.dr, native.dr, windir)
}


## awk

## To use awk to convert a Windows file to Unix, at the Unix prompt, enter:
## awk '{ sub("\r$", ""); print }' winfile.txt > unixfile.txt

## To convert a Unix file to Windows using awk, at the command line, enter:
## awk 'sub("$", "\r")' unixfile.txt > winfile.txt

## On some systems, the version of awk may be old and not include the function sub. If so, try the same command, but with gawk or nawk replacing awk.
## Perl

## To convert a Windows text file to a Unix text file using Perl, at the Unix shell prompt, enter:
## perl -p -e 's/\r$//' < winfile.txt > unixfile.txt

## To convert from a Unix text file to a Windows text file with Perl, at the Unix shell prompt, enter:
## perl -p -e 's/\n/\r\n/' < unixfile.txt > winfile.txt

## unix2dos <- function(unix) {
##   ## this function somehow does not work as expected.
##   ## why? typing the same command in a shell works though.
##   tmp <- tempfile("tmp")
##   on.exit(unlink(tmp))
##   command <- paste("perl -p -e 's/\n/\r\n/' <", unix, " > ", tmp)
##   foo <- system(command)
##   val <- file.copy(tmp, unix, overwrite=TRUE)
##   val
## }

trLbr <- function(unix) {
  lines <- readLines(unix)
  #newlines <- sub('$', '\r', lines)
  #writeLines(newlines, unix)
  writeLines(lines, unix, sep="\r\n")
}

filePathAsAbsolute <- function (x) {
#  if (!file.exists(epath <- path.expand(x))) 
#    stop(gettextf("file '%s' does not exist", x), domain = NA)
  epath <- path.expand(x)
  cwd <- getwd()
  on.exit(setwd(cwd))
  if (tools:::file_test("-d", epath)) {
    setwd(epath)
    getwd()
  }
  else {
    setwd(dirname(epath))
    file.path(getwd(), basename(epath))
  }
}