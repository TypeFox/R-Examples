
#-----------------------------------------------------------------------------
# DateSet

DateSet <- function(name, dates, envir=.Dating, overwrite=FALSE) {
  stopifnot(inherits(name, "character"))
  stopifnot(length(name)==1)
  stopifnot(inherits(dates, c("Date", "POSIXt")))
  stopifnot(!exists(name, envir=envir) | overwrite)
  dating <- new_Dating(c("DateSet", name))
  attr(dating, "environment") <- envir
  assign(paste(name, "Dates", sep="."), dates, envir=envir)
  assign(name, dating, envir=envir)
}

.DDdates <- function(dating) {
  environment <- attr(dating, "environment")
  if(is.null(environment)) environment <- .Dating
  get(paste(dating[length(dating)], "Dates", sep="."), 
    envir=environment)
}

Dbelong.DateSet <- function(dte, dating) {
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  stopifnot(inherits(dating, "DateSet"))
  as.Date(dte) %in% .DDdates(dating)
}

Dseq.DateSet <- function(from, to, dating, len) {
  stopifnot(inherits(dating, "DateSet"))
  stopifnot(inherits(from, c("Date", "POSIXt")))
  alld <- .DDdates(dating)
  posF <- Position(function(x) { x >= as.Date(from) }, alld)
  if(is.na(posF)) return(numeric())
  if(missing(to)) {
    alld[posF:(posF+len-1)]
  } else {
    stopifnot(inherits(to, c("Date", "POSIXt")))
    stopifnot(as.Date(from) <= as.Date(to))
    posT <- Position(function(x) { x > as.Date(to) }, alld)
    if(is.na(posT)) alld[posF:length(alld)]
    else if(posT==1) numeric()
    else alld[posF:(posT-1)]
  }
}

Dsucc.DateSet <- function(dte, dating, num=1) {
  stopifnot(inherits(dating, "DateSet"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  stopifnot(inherits(num, c("numeric", "integer")))
  stopifnot(is.finite(num))
  stopifnot(round(num)==num)
  dte <- as.Date(dte)
  set <- c()
  for(i in 1:length(dte)) {
    nw <- if(num>0) {
      alld <- .DDdates(dating)
      pos <- Position(function(x) { x > as.Date(dte[i]) }, alld)
      if(is.na(pos)) NA
      else if(pos==1) NA
      else if(pos+num-1>length(alld)) NA
      else alld[pos+num-1]
    } else if(num<0) {
      alld <- .DDdates(dating)
      pos <- Position(function(x) { x >= as.Date(dte[i]) }, alld)
      if(is.na(pos)) NA
      else if(pos+num<1) NA
      else alld[pos+num]
    } else {
      Dround(dte[i], dating)
    }
    set <- append(set, nw)
  }
  set
}

Dfloor.DateSet <- function(dte, dating) {
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  stopifnot(inherits(dating, "DateSet"))
  dte <- as.Date(dte)
  alld <- .DDdates(dating)
  set <- c()
  for(k in 1:length(dte)) {
    pos <- Position(function(x) { x > dte[k] }, alld)
    nw <- if(is.na(pos)) alld[length(alld)]
    else if(pos==1) as.Date(NA)
    else alld[pos-1]
    set <- append(set, nw)
  }
  set
}

Dceiling.DateSet <- function(dte, dating) {
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  stopifnot(inherits(dating, "DateSet"))
  dte <- as.Date(dte)
  alld <- .DDdates(dating)
  set <- c()
  for(k in 1:length(dte)) {
    pos <- Position(function(x) { x >= dte[k] }, alld)
    nw <- if(is.na(pos)) as.Date(NA)
    else alld[pos]
    set <- append(set, nw)
  }
  set
}

Dround.DateSet <- function(dte, dating) {
  stopifnot(inherits(dating, "DateSet"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  alld <- .DDdates(dating)
  set <- c()
  for(k in 1:length(dte)) {
    pos <- Position(function(x) { x >= dte[k] }, alld)
    nw <- if(is.na(pos)) alld[length(alld)] #Dfloor
    else if(pos==1) alld[1] #Dceiling
    else if(alld[pos]==dte[k]) alld[pos]
    else if(Ddiff(alld[pos-1],dte[k],Daily)>Ddiff(dte[k],alld[pos],Daily)) alld[pos] 
    else alld[pos-1]
    set <- append(set, nw)
  }
  set
}

Ddiff.DateSet <- function(dte1, dte2, dating) {
  stopifnot(inherits(dating, "DateSet"))
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dte1 <- as.Date(dte1)
  dte2 <- as.Date(dte2)
  alld <- .DDdates(dating)
  pos1 <- Position(function(x) { x > as.Date(dte1) }, alld)
  if(is.na(pos1)) return(NA)
  if(pos1==1) return(NA)
  pos2 <- Position(function(x) { x > as.Date(dte2) }, alld)
  if(is.na(pos2)) return(NA)
  if(pos2==1) return(NA)
  pos2-pos1
}

#-----------------------------------------------------------------------------
