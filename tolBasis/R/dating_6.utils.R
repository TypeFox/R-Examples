
as.date.list <- function(dates) {
  if(inherits(dates, "Date")) return(dates)
  else if(inherits(dates, "POSIXt")) {
    if(all(Dbelong(dates, Daily))) return(Dround(dates, Daily))
    else stop("Dates should not have fractional part.")
  } else stop("Dates should inherit classes Date or POSIXt.")
}

is.date.sequence <- function(date.sequence) {
  stopifnot(inherits(date.sequence, "Date"))
  all(convolve(as.numeric(date.sequence), c(-1,1), type="f")>0)
}

.Dcheck <- function(date.sequence, dating) {
  mindate <- min(date.sequence)
  if(!Dbelong(mindate, dating)) return(FALSE)
  maxdate <- max(date.sequence)
  if(!Dbelong(maxdate, dating)) return(FALSE)
  len <- Ddiff(mindate, maxdate, dating)
  if(len+1==length(date.sequence)) {
    dates2 <- Dseq(mindate, maxdate, dating)
    setequal(date.sequence, dates2)
  } else FALSE
}

Dcheck <- function(date.sequence, dating) {
  date.sequence <- as.date.list(date.sequence) # stop if is not a date list
  stopifnot(is.date.sequence(date.sequence))
  stopifnot(inherits(dating, "Dating"))
  .Dcheck(date.sequence, dating)
}

Dfind <- function(date.sequence) {
  date.sequence <- as.date.list(date.sequence) # stop if is not a date list
  stopifnot(is.date.sequence(date.sequence))
  datings <- Datings()
  for(i in 1:length(datings))
    if(.Dcheck(date.sequence, datings[[i]])) return(datings[[i]])
  return(NULL)
}

# Datings(): Obtain a list with the defined datings
Datings <- function() {
  # predefined datings
  dd1 <- list(Yearly, HalfYearly, Quarterly, Monthly, Mondays, Tuesdays, 
    Wednesdays, Thursdays, Fridays, Saturdays, Sundays, Daily)
  # user-defined datesets
  dlist <- ls(envir=.Dating)
  objects <- lapply(dlist, function(x) { eval(parse(text=x), envir=.Dating) })
  check <- lapply(objects, function(x) { inherits(x, "Dating")})
  dd2 <- objects[unlist(check)]
  # defined datings
  c(dd1, dd2)
}
