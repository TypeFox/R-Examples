
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

.is.numeric_wholenumber_not0 <- function(x, tol = .Machine$double.eps^0.5) {
  if(is.numeric(x)) {
    if(.is.wholenumber(x, tol)) {
      round(x)!=0
    } else FALSE
  } else FALSE
}

#-----------------------------------------------------------------------------

Serie <- function(data, dating, begin) {
  stopifnot(inherits(data, c("numeric", "integer")))
  stopifnot(length(data)>0)
  stopifnot(inherits(dating, "Dating"))
  stopifnot(inherits(begin, c("Date", "POSIXt")))
  stopifnot(Dbelong(begin, dating))
  class(data) <- c("Serie")
  attr(data, "dating") <- dating
  attr(data, "begin") <- as.Date(begin)
  data
}

#-----------------------------------------------------------------------------

Sdating <- function(s) UseMethod("Sdating")

Sdating.Serie <- function(s) attr(s, "dating")

#-----------------------------------------------------------------------------

print.Serie <- function(x, ..., limit, mode) {
  default <- missing(limit) & missing(mode)
  if(missing(limit)) limit <- getOption("max.print.Serie")
  if(missing(mode)) mode <- getOption("mode.print.Serie")
  stopifnot(is.na(limit) | (limit>=1 & round(limit)==limit))
  stopifnot(is.numeric(mode))
  limited <- if(is.na(limit)) FALSE else (length(x)>limit)
  range1 <- if(!limited) 1:length(x)
    else if(mode<0) 1:limit
    else if(mode>0) numeric()
    else 1:floor(limit/2)
  range2 <- if(!limited) numeric()
    else if(mode<0) numeric()
    else if(mode>0) (length(x)-limit+1):length(x)
    else (length(x)-ceiling(limit/2)+1):length(x)
  dtes <- Dseq(attr(x, "begin"), dating=Sdating(x), len=length(x))
  if(length(range1)) for(i in range1) {
    cat(as.character(dtes[i]), x[i], "\n")
  } 
  if(limited) cat(paste0("... (", length(x)-limit, " more rows)\n"))
  if(length(range2)) for(i in range2) {
    cat(as.character(dtes[i]), x[i], "\n")
  } 
  if(limited & default) 
    cat("[ see options 'limit' and 'mode' at help(print.Serie) ]\n")
  invisible(x) 
}

#-----------------------------------------------------------------------------

Sfirst <- function(s) UseMethod("Sfirst")
Sfirst.Serie <- function(s) {
  stopifnot(inherits(s, "Serie"))
  attr(s, "begin")
}

Slast <- function(s) UseMethod("Slast")
Slast.Serie <- function(s) {
  stopifnot(inherits(s, "Serie"))
  Dsucc(attr(s, "begin"), Sdating(s), length(s)-1)
}

Ssub <- function(s, from=NA, to=NA) UseMethod("Ssub")
Ssub.Serie <- function(s, from=NA, to=NA) {
  stopifnot(inherits(s, "Serie"))
  stopifnot(is.na(from) | .is.numeric_wholenumber_not0(from) | inherits(from, c("Date", "POSIXt")))
  stopifnot(is.na(to) | .is.numeric_wholenumber_not0(to) | inherits(to, c("Date", "POSIXt")))
  ini1 <- Sfirst(s)
  fin1 <- Slast(s)
  ini2 <- if(is.na(from)) ini1 
    else if(.is.numeric_wholenumber_not0(from)) { 
      if(from>0) Dsucc(ini1, Sdating(s), from-1)
      else Dsucc(fin1, Sdating(s), from+1)
    } else from
  fin2 <- if(is.na(to)) fin1 
    else if(.is.numeric_wholenumber_not0(to)) {
      if(to>0) Dsucc(ini1, Sdating(s), to-1)
      else Dsucc(fin1, Sdating(s), to+1)
    } else to
  if(Ddiff(ini2,fin2,Sdating(s))<0)
    stop("The argument 'from' should be earlier than the argument 'to'.")
  i1 <- Ddiff(ini1,ini2,Sdating(s)) + 1
  i2 <- Ddiff(ini1,fin2,Sdating(s)) + 1
  c0 <- if(ini2<ini1) rep(NA,Ddiff(ini2,ini1,Sdating(s))) else numeric()
  c1 <- s[max(1,i1):min(i2,length(s))]
  c2 <- if(fin1<fin2) rep(NA,Ddiff(fin1,fin2,Sdating(s))) else numeric()
  Serie(c(c0,c1,c2), Sdating(s), ini2)
}

Sdates <- function(s) UseMethod("Sdates")
Sdates.Serie <- function(s) {
  stopifnot(inherits(s, "Serie"))
  Dseq(attr(s, "begin"), dating=attr(s, "dating"), len=length(s))
}

#-----------------------------------------------------------------------------

`+.Serie` <- function(s1, s2) {
  stopifnot(inherits(s1, c("Serie", "numeric", "integer")))
  stopifnot(inherits(s2, c("Serie", "numeric", "integer")))
  if(inherits(s1, c("numeric", "integer"))) {
    ss <- s2
    ss[1:length(ss)] <- ss[1:length(ss)] + s1
    ss
  } else if(inherits(s2, c("numeric", "integer"))) {
    ss <- s1
    ss[1:length(ss)] <- ss[1:length(ss)] + s2
    ss
  } else {
    d1 <- Sdating(s1)
    d2 <- Sdating(s2)
    if(identical(d1,d2)) {
      ini <- max(Sfirst(s1), Sfirst(s2))
      fin <- min(Slast(s1), Slast(s2))
      ss <- Serie(as.numeric(Ssub(s1, ini, fin))+as.numeric(Ssub(s2, ini, fin)), d1, ini)
    } else {
      # different datings
      stop()
    }
  }
}

`-.Serie` <- function(s1, s2) {
  stopifnot(inherits(s1, c("Serie", "numeric", "integer")))
  stopifnot(inherits(s2, c("Serie", "numeric", "integer")))
  s1 + (s2*(-1))
}

`*.Serie` <- function(s1, s2) {
  stopifnot(inherits(s1, c("Serie", "numeric", "integer")))
  stopifnot(inherits(s2, c("Serie", "numeric", "integer")))
  if(inherits(s1, c("numeric", "integer"))) {
    ss <- s2
    ss[1:length(ss)] <- ss[1:length(ss)] * s1
    ss
  } else if(inherits(s2, c("numeric", "integer"))) {
    ss <- s1
    ss[1:length(ss)] <- ss[1:length(ss)] * s2
    ss
  } else {
    d1 <- Sdating(s1)
    d2 <- Sdating(s2)
    if(identical(d1,d2)) {
      ini <- max(Sfirst(s1), Sfirst(s2))
      fin <- min(Slast(s1), Slast(s2))
      ss <- Serie(as.numeric(Ssub(s1, ini, fin))*as.numeric(Ssub(s2, ini, fin)), d1, ini)
    } else {
      # different datings
      stop()
    }
  }
}

`/.Serie` <- function(s1, s2) {
  stopifnot(inherits(s1, c("Serie", "numeric", "integer")))
  stopifnot(inherits(s2, c("Serie", "numeric", "integer")))
  s1 * (s2^(-1))
}

`^.Serie` <- function(s, n) {
  stopifnot(inherits(s, "Serie"))
  stopifnot(inherits(n, c("numeric", "integer")))
  Serie(as.numeric(s)^n, Sdating(s), Sfirst(s))
}

`[.Serie` <- function(s, index) {
  stopifnot(inherits(s, "Serie"))
  if(missing(index)) return(s)
  snum <- as.numeric(s)
  if(inherits(index, c("Date", "POSIXt"))) {
    dates <- as.Date(index)
    b <- attr(s, "begin")
    dating <- attr(s, "dating")
    sapply(dates, function(d) { 
      if(Dbelong(d, dating)) {
        i <- Ddiff(b, d, dating) + 1
        if(i<1 | i > length(s)) NA else snum[i]
      } else NA
    })
  } else snum[index]
}

#`%at%` <- function(s1, s2) { Serie(s1[Sdates(s2)], Sdating(s2), Sfirst(s2)) }
#lm(sx ~ B %:% sx %at% sx)

#-----------------------------------------------------------------------------

as.Serie <- function(x, ...) UseMethod("as.Serie")

as.Serie.ts <- function(x, ...) {
  stopifnot(inherits(x, "ts"))
  sinfo <- attr(x, "tsp")
  if(sinfo[3]==1 & .is.wholenumber(sinfo[1])) {
    begin <- as.Date(ISOdate(sinfo[1],1,1))
    dating <- Yearly
  } else if(sinfo[3]==4 & .is.wholenumber(sinfo[1]*sinfo[3])) {
    begin <- as.Date(ISOdate(floor(sinfo[1]), 1+3*round(sinfo[1]*sinfo[3])%%sinfo[3], 1))
    dating <- Quarterly
  } else if(sinfo[3]==12 & .is.wholenumber(sinfo[1]*sinfo[3])) {
    begin <- as.Date(ISOdate(floor(sinfo[1]), 1+round(sinfo[1]*sinfo[3])%%sinfo[3], 1))
    dating <- Monthly
  } else stop()
  if(inherits(x, "mts")) {
    listS <- lapply(1:ncol(x), function(i) { Serie(as.numeric(x[,i]), dating, begin) })
    names(listS) <- colnames(x)
    listS     
  } else {
    Serie(as.numeric(x), dating, begin)
  }
}

as.ts.Serie <- function(x, ...) {
  freq <- Ddiff(as.Date(ISOdate(1970,1,1)), as.Date(ISOdate(2070,1,1)), Sdating(x)) / 100
  dif <- Ddiff(as.Date(ISOdate(2000,1,1)), Sfirst(x), Sdating(x))
  ts(as.numeric(x), start = 2000 + dif/freq, frequency = freq)
}

as.Serie.xts <- function(x, ..., dating) {
  stopifnot(requireNamespace("zoo", quietly=TRUE))
  stopifnot(inherits(x, "xts"))
  if(missing(dating)) {
    dating <- Dfind(zoo::index(x))
    if(is.null(dating)) stop("Cannot find a compatible dating.")
  } else {
    stopifnot(inherits(x, "Dating"))
    if(!all(Dbelong(zoo::index(x), dating))) stop("The indicated dating is not compatible.")
  }
  begin <- as.Date(min(zoo::index(x)))
  if(ncol(x)==1) Serie(as.numeric(x), dating, begin)
  else {
    listS <- lapply(1:ncol(x), function(i) { Serie(as.numeric(x[,i]), dating, begin) })
    names(listS) <- colnames(x)
    listS
  }
}

as.xts.Serie <- function(x, ...) {
  stopifnot(requireNamespace("xts", quietly=TRUE))
  xts::xts(as.matrix(x), order.by=Sdates(x))
}

#-----------------------------------------------------------------------------

`%:%` <- function(p, s) {
  stopifnot(inherits(p, "Polyn"))
  stopifnot(inherits(s, "Serie"))
  plen <- length(p)
  slen <- length(s)
  if(plen>slen) stop() else {
    dn <- rep(0, slen-plen+1)
    for(i in 1:plen) {
      dn <- dn + as.numeric(s)[i:(slen-plen+i)] * p[i]
    }
    Serie(dn, Sdating(s), Dsucc(Sfirst(s), Sdating(s), attr(p, "base")+plen-1))
  }
}

#-----------------------------------------------------------------------------
