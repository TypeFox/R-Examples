
#-----------------------------------------------------------------------------
# NMonthly

Quarterly <- new_Dating(c("NMonthly", "Quarterly"))
HalfYearly <- new_Dating(c("NMonthly", "HalfYearly"))

.DDmonthfrequency <- function(x) UseMethod(".DDmonthfrequency")
.DDmonthfrequency.Quarterly <- function(s) { 3 }
.DDmonthfrequency.HalfYearly <- function(s) { 6 }

Dbelong.NMonthly <- function(dte, dating) {
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  stopifnot(inherits(dating, "NMonthly"))
  freq <- .DDmonthfrequency(dating) 
  ref <- as.Date(ISOdate(2000,1,1))
  (Ddiff(ref, as.Date(dte), Monthly) %% freq) == 0
}

Dseq.NMonthly <- function(from, to, dating, len) {
  stopifnot(inherits(dating, "NMonthly"))
  stopifnot(inherits(from, c("Date", "POSIXt")))
  from <- as.Date(from)
  freq <- .DDmonthfrequency(dating)
  if(missing(to)) {
    seq(as.Date(Dceiling(from, dating)), length.out=len, by=paste(freq, "months"))
  } else {
    stopifnot(inherits(to, c("Date", "POSIXt")))
    to <- as.Date(to)
    stopifnot(from <= to)
    seq(as.Date(Dceiling(from, dating)), as.Date(Dfloor(to, dating)), paste(freq, "months"))
  }
}

Dsucc.NMonthly <- function(dte, dating, num=1) {
  stopifnot(inherits(dating, "NMonthly"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  freq <- .DDmonthfrequency(dating)
  Dsucc(dte, Monthly, num*freq)
}

Dfloor.NMonthly <- function(dte, dating) {
  stopifnot(inherits(dating, "NMonthly"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  freq <- .DDmonthfrequency(dating)
  dteM <- Dfloor(dte, Monthly)
  prev <- (Ddiff(as.Date(ISOdate(2000,1,1)), dteM, Monthly) %% freq) 
  # ifelse(prev>0, Dsucc(dteM, Monthly, -prev), dteM)
  set <- c()
  for(i in 1:length(prev)) {
    nw <- if(prev[i]>0) Dsucc(dteM[i], Monthly, -prev[i]) else dteM[i]
    set <- append(set, nw)
  }
  set
}

Dceiling.NMonthly <- function(dte, dating) {
  stopifnot(inherits(dating, "NMonthly"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  freq <- .DDmonthfrequency(dating)
  dteM <- Dceiling(dte, Monthly)
  post <- (Ddiff(as.Date(ISOdate(2000,1,1)), dteM, Monthly) %% freq) 
  #ifelse(post>0, Dsucc(dteM, Monthly, freq-post), dteM)
  set <- c()
  for(i in 1:length(post)) {
    nw <- if(post[i]>0) Dsucc(dteM[i], Monthly, freq-post[i]) else dteM[i]
    set <- append(set, nw)
  }
  set
}

Dround.NMonthly <- function(dte, dating) {
  stopifnot(inherits(dating, "NMonthly"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  dteF <- Dfloor.NMonthly(dte, dating)
  dteC <- Dceiling.NMonthly(dte, dating)
  ifelse(Ddiff(dteF,dte,Daily)>Ddiff(dte,dteC,Daily), dteC, dteF)
}

Ddiff.NMonthly <- function(dte1, dte2, dating) {
  stopifnot(inherits(dating, "NMonthly"))
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dte1 <- as.Date(dte1)
  dte2 <- as.Date(dte2)
  freq <- .DDmonthfrequency(dating)
  Ddiff(Dfloor(dte1, dating), Dfloor(dte2, dating), Monthly) / freq
}

#-----------------------------------------------------------------------------
