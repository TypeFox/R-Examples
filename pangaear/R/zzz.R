baseoai <- function() "http://ws.pangaea.de/oai/"
base <- function() 'https://doi.pangaea.de/'
sbase <- function() "https://www.pangaea.de/search"

pgc <- function(x) Filter(Negate(is.null), x)

pluck <- function(x, name, type) {
  if (missing(type)) {
    lapply(x, "[[", name)
  } else {
    vapply(x, "[[", name, FUN.VALUE = type)
  }
}

check <- function(x){
  if(is.character(x)){
    if( grepl("does not exist|unknown", x))
      stop(x, call. = FALSE)
  }
}

read_csv <- function(x){
  lns <- readLines(x, n = 300)
  ln_no <- grep("\\*/", lns)
  tmp <- read.csv(x, header = FALSE, sep = "\t", skip = ln_no+1, stringsAsFactors=FALSE)
  nn <- strsplit(lns[ln_no+1], "\t")[[1]]
  setNames(tmp, nn)
}

ifn <- function(x) if(is.null(x)) NA else x

strextract <- function(str, pattern) regmatches(str, regexpr(pattern, str))

capwords <- function(s, strict = FALSE, onlyfirst = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)), {
    s <- substring(s,2); if(strict) tolower(s) else s
  }, sep = "", collapse = " " )
  if(!onlyfirst){
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  } else {
    sapply(s, function(x)
      paste(toupper(substring(x,1,1)),
            tolower(substring(x,2)),
            sep="", collapse=" "), USE.NAMES=F)
  }
}
