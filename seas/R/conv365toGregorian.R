conv365toGregorian <-
function(x) {
  sc <- seas.df.check(x)  # make sure it has a date column
  a <- list()  # make a copy of the column attributes
  for (n in names(x)) {
    tmp <- attributes(x[[n]])  # strip out unneeded attributes
    srp <- names(tmp) %in% c("dim", "dimnames", "names")
    a[[n]] <- if (any(srp))
      tmp[!srp]
    else
      tmp
  }
  if (any(is.na(x)))
    warning("'x' contain 'NA' values which may be repeated unexpectedly!")
  date.range <- range(x$date)
  if (as.integer(format(date.range[2], "%d")) == 30)
    date.range[2] <- date.range[2] + 1
  Gdates <- data.frame(date=seq(date.range[1], date.range[2], by="day"))
  comb <- merge(x, Gdates, by="date", all.y=TRUE)
  nm <- names(comb)
  nm <- nm[nm != "date"]
  miss <- which(apply(comb[,nm, drop=FALSE], 1, function(x) all(is.na(x))))
  if (any(diff(miss) <= 1))
    stop("More than one consecutive day(s) missing")
  comb[miss,] <- comb[miss - 1,]  # the day before
  for (n in names(comb))  # copy attributes back
    attributes(comb[[n]]) <- a[[n]]
  comb
}
