# moving average as described in Tallaksen and van Lanen (2004)
# where the n past values are averaged
ma <- function(x, n, sides = 1)  {
  y <- filter(x, filter = rep(x = 1/n, times = n), sides = sides)


  # filter() returns a ts-object
  return(as.numeric(y))
}

.check_xts <- function(x) {
  # check if a regular time series is provided
  dt <- diff(time(x))
  dt <- dt[!duplicated(dt)]    # unique() drops unit
  if(length(dt) != 1) stop("Only regular time series are supported.")

  dt <- as.numeric(dt, units = "secs")
  xtsAttributes(x)[["deltat"]] <- dt

  # check if a unit is provided
  unit <- xtsAttributes(x)[["unit"]]
  if(is.null(unit) || unit == "" || is.na(unit)) {
    stop("no unit found in attributes")
  } else {
    # if so, parse volume und time
    names(unit) <- "flow"
    xtsAttributes(x)[["unit.parsed"]] <- .split_unit(unit)
  }

  # set colnames
  if(ncol(x) == 1) colnames(x) <- "discharge"

  return(x)
}

.split_unit <- function(x) {
  y <- strsplit(gsub("\\^.", "", x), "/")[[1]]
  y <- gsub("\u00B3", "", y)
  dict <- c("s" = "secs", "m" = "mins", "h" = "hours", "d" = "days")

  units <- c(volume = y[1], time = unname(dict[y[2]]))
  if(is.na(units["time"])) stop("unknown time unit ", sQuote(y[2]))

  return(units)
}


.conv_factor <- function(from, to, dimension = c("time", "volume")) {
  dict <- list(time = c("days" = 86400, "hours" = 3600, "mins" = 60, "secs" = 1),
               volume = c("m" = 1, "l" = 1e-3, "cm" = 1e-6))

  dimension <- match.arg(dimension)
  from <- match.arg(from, names(dict[[dimension]]), several.ok = F)
  to <- match.arg(to, names(dict[[dimension]]), several.ok = F)

  x <- unname(dict[[dimension]][to]/dict[[dimension]][from])
  return(x)
}


as.xts.lfobj <- function(x, ...) {
  lfcheck(x)

  time <- with(x, as.Date(paste(year, month, day, sep = "-")))
  y <- xts(x[, "flow"], order.by = time)

  att <- attr(x, "lfobj")
  #att[["location"]] <- att[["station"]]
  missing <- setdiff(c("river", "station", "unit", "institution"), names(att))
  att[missing] <- ""
  xtsAttributes(y) <- att

  # correct colname is set by .check_xts
  y <- .check_xts(y)

  return(y)
}


# classify values due to their neighbours
group <- function(x, new.group.na = TRUE) {
  inc <- diff(as.numeric(x))
  if (new.group.na) inc[is.na(inc)] <- Inf

  grp <- c(0, cumsum(inc != 0))

  if(grp[1] == 0) grp <- grp + 1

  return(factor(grp))
}




period <- function(x, varying) {
  x <- as.xts(x)

  if (is.character(varying)){
    if(varying == "constant") return(rep_len(1, nrow(x)))

    # periodically varying threshold
    f <- c("daily" = "%j", "weekly" = "%V", "monthly" = "%m")
    varying <- match.arg(arg = varying, choices = names(f), several.ok = FALSE)

    period <- as.numeric(format(time(x), format = f[varying]))
  } else {
    # is.Date(varying): seasonal threshold
    varying <- as.Date(varying)
    day <- as.numeric(format(time(x), "%j"))
    ep <- sort(as.numeric(format(varying, format = "%j")))

    period <- rep(ep[1], nrow(x))
    for(i in rev(ep)) period[day < i] <- i
  }

  return(period)
}


vary_threshold <- function(x, varying = "constant",
                           fun = function(x)
                             quantile(x, probs = 0.05, na.rm = TRUE),
                           ...) {
  x <- as.xts(x)

  zz <- x
  g <- period(x, varying = varying)
  split(zz, g) <- lapply(split(coredata(x), g), FUN = fun, ...)

  colnames(zz) <- "threshold"

  return(zz)
}



strsplit_date <- function(x, prefix = "") {
  time <- as.Date(x)

  format <- c(day = "%d", month = "%m", year = "%Y")
  y <- lapply(format, function(f) as.numeric(format(time, format = f)))
  y <- do.call(cbind, y)

  colnames(y) <- paste0(prefix, colnames(y))

  return(y)
}

# hack, because in all.equal() the user can't enforce the interpretation of
# argument tolerance as absolute differences
expect_equal2 <- function(object, expected, tolerance = 1e-10, ...) {
  testthat::expect_true(all(abs(object - expected) < tolerance), ...)
}


.char2html <- function(x, dict = c("^2" = "sup2",  "^3"= "sup3",
                                   "\u00df"= "szlig", # scharfes S
                                   "\u00e4" = "auml", "\u00c4" = "Auml", # ae AE
                                   "\u00f6" = "ouml", "\u00d6" = "Ouml", # oe OE
                                   "\u00fc" = "uuml", "\u00dc" = "Uuml"  # ue UE
)) {
  tbl <- paste0("&", dict, ";")
  names(tbl) <- names(dict)

  for(i in seq_along(tbl)){
    x <- gsub(names(tbl[i]), tbl[i], x, fixed = T)
  }

  return(x)
}
