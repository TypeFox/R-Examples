
#-----------------------------------------------------------------------------
# Date

Date <- function(year, month=1, day=1) {
  #as.Date(ymd(year*10000+month*100+day))
  as.Date(strptime(paste(year, month, day, sep="-"), "%Y-%m-%d", tz="GMT"))
}

#-----------------------------------------------------------------------------
# Dating

new_Dating <- function(name) {
  stopifnot(inherits(name, "character"))
  dating <- name[length(name)]
  class(dating) <- c("Dating", name)
  dating
}

Daily <- new_Dating("Daily")

print.Dating <- function(x,...) {
  print(as.character(x))
  invisible(x)
}

#-----------------------------------------------------------------------------
# generic

Dbelong <- function(dte, dating=Daily) UseMethod("Dbelong", dating)
Dseq <- function(from, to, dating, len) UseMethod("Dseq", dating)

Dsucc <- function(dte, dating=Daily, num=1) UseMethod("Dsucc", dating)
Dfloor <- function(dte, dating=Daily) UseMethod("Dfloor", dating)
Dceiling <- function(dte, dating=Daily) UseMethod("Dceiling", dating)
Dround <- function(dte, dating=Daily) UseMethod("Dround", dating)
Ddiff <- function(dte1, dte2, dating=Daily) UseMethod("Ddiff", dating)

#-----------------------------------------------------------------------------
