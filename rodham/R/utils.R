raw2clean <- function(x){
  temp <- splitstackshape::cSplit(x, "to", "; ", "long")
  clean <- splitstackshape::cSplit(temp, "from", "; ", "long")
  return(clean)
}

checkRelease <- function(x){
  x <- tolower(x) # norm
  if (x %in% db$name && length(x) == 1) {
    release <- db[db$name == x, "uri"]
    return(release)
  } else if(!x %in% db$name && length(x) == 1) {
    msg <- paste0(db$name, collapse = ", ")
    stop(paste0("Wrong release, can be one of: ", msg), call. = FALSE)
  } else if (length(x) > 1){
    stop("Can only pass one release", call. = FALSE)
  }
}

OStoURI <- function(x){
  if(x == "Windows"){
    uri <- "ftp://ftp.foolabs.com/pub/xpdf/xpdfbin-win-3.04.zip"
    ext <- ".zip"
  } else if (x == "Linux"){
    uri <- "ftp://ftp.foolabs.com/pub/xpdf/xpdfbin-linux-3.04.tar.gz"
    ext <- ".tar.gz"
  } else if(x == "Darwin") {
    uri <- "ftp://ftp.foolabs.com/pub/xpdf/xpdfbin-mac-3.04.tar.gz"
    ext <- ".tar.gz"
  } else {
    stop(paste0("Unrecognised operating system, download manually from \n",
                "http://www.foolabs.com/xpdf/download.html"), call. = FALSE)
  }
  lst <- list(uri = uri, ext = ext)
  return(lst)
}
