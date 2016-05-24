"seas.sum" <-
function(x, var, width = 11, start.day = 1,
         prime, a.cut = 0.3, na.cut = 0.2) {
  orig <- as.character(substitute(x))[[1]]
  if (missing(var))  # try and find something useful
    var <- c("precip", "rain", "snow", "leak", "evap", "ezwat", "et", "runoff", "air", "soil")
  var <- names(x)[names(x) %in% var]
  if (missing(prime))  # use the first variable as 'prime'
    prime <- ifelse("precip" %in% var, "precip", var[1])
  sc <- seas.df.check(x, orig, c(prime, var))
  if (is.na(a.cut) || a.cut <= 0)
    a.cut = FALSE
  x$fact <- mkseas(x, width, start.day)
  x$ann <- mkann(x, start.day)
  bins <- levels(x$fact)
  num <- length(bins)
  years <- levels(x$ann)
  ann <- data.frame(year=years, active=NA, days=NA, na=NA)
  seas <- array(dim=c(length(years), num, length(var)),
                dimnames=list(years, bins, var))
  days <- array(dim=c(length(years), num),
                dimnames=list(years, bins))
  na <- days
  if (is.na(a.cut) || !a.cut) {
    a.cut <- FALSE
  } else {
    active <- seas  # copy array
    is.active <- function(test) {  # test to count the number of active days
      tot <- numeric(length(test))
      tot[is.na(test)] <- NA  # keep NAs
      tot[test > a.cut] <- 1  # find only days where the prime > a.cut
      na.rm <- ifelse(sum(is.na(tot)) / length(tot) < na.cut[2], TRUE, FALSE)
      return(sum(tot, na.rm=na.rm))
    }
  }
  if (length(na.cut) == 1)
    na.cut <- rep(na.cut, 2)
  sum.is.num <- function(x)sum(is.finite(x))
  ann$days <- attr(x$ann, "year.lengths")
  if (a.cut)
    ann$active <- tapply(x[,prime], x$ann, is.active)
  else
    ann$active <- NULL
  ann$na <- tapply(x[,prime], x$ann, sum.is.num)
  for (p in var)
    ann[,p] <- tapply(x[,p], x$ann, sum, na.rm=TRUE)
  td <- function(y) table(mkseas(width=width, year=y, calendar=sc$calendar))
  days[,] <- t(sapply(ann$days, td))
  for (y in 1:length(years)) {  # using integer indexes
    s <- x[x$ann == years[y],, drop=FALSE]
    if (nrow(s) > 0) {
      na[y,] <- tapply(s[,prime], s$fact, sum.is.num)
      for (p in var) {
        seas[y,,p] <- tapply(s[,p], s$fact, sum, na.rm=TRUE)
        if (a.cut)
          active[y,,p] <- tapply(s[,p], s$fact, is.active)
      }
    }
  }
  ann$na[is.na(ann$na)] <- 0
  ann$na <- ann$days - ann$na
  na[is.na(na)] <- 0
  na <- days - na
  ann.na <- ann$na / ann$days > na.cut[1]
  seas.na <- na/days > na.cut[2]
  ann[ann.na, var] <- NA
  seas[,,var][seas.na] <- NA
  if (a.cut) {
    ann[ann.na, "active"] <- NA
    active[,,var][seas.na] <- NA
  }
  l <- list(ann=ann, seas=seas, active=NA, days=days, na=na)
  if (a.cut)
    l$active <- active
  else
    l$active <- NULL
  l$start.day <- start.day
  l$years <- years
  l$var <- var
  l$units <- list()
  l$long.name <- list()
  for (v in var) {
    l$units[[v]] <- attr(x[[v]], "units")
    l$long.name[[v]] <- attr(x[[v]], "long.name")
    if (is.null(l$long.name[[v]]))
      l$long.name[[v]] <- v  # use var name
  }
  l$prime <- prime
  l$width <- width
  l$bins <- bins
  l$bin.lengths <- attr(x$fact, "bin.lengths")
  l$year.range <- attr(x$ann, "year.range")
  l$na.cut <- na.cut
  l$a.cut <- a.cut
  l$id <- sc$id
  l$name <- sc$name
  attr(l, "class") <- "seas.sum"
  return(l)
}
