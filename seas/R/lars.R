"read.lars" <-
function(stfile, year.offset = 0) {
  conv <- data.frame(
    lars=I(c("year", "jday", "rad", "min", "max", "rain", "sun")),
    seas=I(c("year", "yday", "solar", "t_min", "t_max", "precip", "sun")))
  if (is.character(stfile)) {
    ext <- unlist(strsplit(stfile, "\\."))
    ext <- tolower(ext[length(ext)])
    if (ext != "st")
      warning("This 'stfile' normally has a '.st' extension")
    ld <- dirname(stfile)
  } else {
    ld <- "."
  }
  st <- readLines(stfile, warn=FALSE)
  datfile <- st[match("[WEATHER FILES]", st) + 1]
  if (is.na(datfile))
    stop("Invalid '[WEATHER FILES]'")
  if (basename(datfile) == datfile)  # make full path
    datfile <- paste(ld, basename(datfile), sep="/")
  format <- tolower(unlist(strsplit(st[match("[FORMAT]", st) + 1],
                                    "[ \t]+", perl=TRUE)))
  m <- data.frame(m=match(format, conv$lars))
  m$i <- 1:nrow(m)
  m <- m[!is.na(m$m),]
  format[m$i] <- conv$seas[m$m]
  x <- read.table(datfile, col.names=format, na.strings="-99")
  x$year <- as.integer(x$year + year.offset)
  x$date <- if ("month" %in% names(x) && !"yday" %in% names(x))
    with(x, as.Date(sprintf("%04i-%02i-%02i", as.integer(year),
                           as.integer(month), as.integer(day))))
  else
    with(x, as.Date(sprintf("%04i-%03i", as.integer(year), as.integer(yday)),
                    "%Y-%j"))
  deg <- iconv("\260", "latin1", "")
  deg <- if (is.na(deg)) "deg" else deg
  degC <- paste(deg, "C", sep="")
  if ("t_max" %in% names(x)) {
    attr(x$t_max, "units") <- degC
    attr(x$t_max, "long.name") <- "daily maximum temperature"
  }
  if ("t_min" %in% names(x)) {
    attr(x$t_min, "units") <- degC
    attr(x$t_min, "long.name") <- "daily minimum temperature"
  }
  if ("precip" %in% names(x)) {
    attr(x$precip, "units") <- "mm"
    attr(x$precip, "long.name") <- "total daily precipitation"
  }
  if ("solar" %in% names(x)) {
    # The LARS-WG documentation says "MJ/day"
    # however, it is more precisely "per m^2" since it is a flux
    sq <- iconv("\262", "latin1", "")
    sq <- if (is.na(sq)) "^2" else sq
    radday <- paste("MJ/(m", sq, "day)", sep="")
    attr(x$solar, "units") <- radday
    attr(x$solar, "long.name") <- "total daily global solar radiation"
  }
  if ("pet" %in% names(x)) {
    attr(x$pet, "units") <- "mm"
    attr(x$pet, "long.name") <- "potential evapotranspiration"
  }
  if ("[SITE]" %in% st) {
    site <- st[match("[SITE]", st) + 1]
    if (site != "")
      attr(x, "name") <- site
  }
  if ("[LAT, LON and ALT]" %in% st) {
    lla <- strsplit(st[match("[LAT, LON and ALT]", st) + 1],
                    "[ \t]+", perl=TRUE)
    lla <- as.numeric(unlist(lla))
    lla <- lla[!is.na(lla)]
    if (length(lla) == 3) {
      attr(x, "latitude") <- lla[1]
      attr(x, "longitude") <- lla[2]
      attr(x, "elevation") <- lla[3]
    }
  }
  x
}

"write.lars" <-
function(x, stfile, datfile, site, lat, lon, alt) {
  conv <- data.frame(
    lars=I(c("year", "jday", "min", "max", "rain", "rad", "sun")),
    seas=I(c("year", "yday", "t_min", "t_max", "precip", "solar", "sun")))
  # pre-proces st file, if selected
  mkst <- TRUE
  if (is.na(stfile) || (is.logical(stfile) && !stfile)) {
    mkst <- FALSE
  } else {
    if (!is.character(stfile))
      stop("'stfile' must be a file pathname to an '.st' file")
    stfilepart <- unlist(strsplit(basename(stfile), "\\."))
    if (tolower(stfilepart[length(stfilepart)]) != "st")
      warning("This 'stfile' normally has a '.st' extension")
    ld <- dirname(stfile)
    if (missing(datfile)) {
      datfile <- sprintf("%s/%s.sr", ld, stfilepart[1])
    }
  }
  if (mkst) {
    checkattr <- function(a, n) {
      if (a %in% names(attributes(x)))
        return(attr(x, a))
      else
        stop(sprintf("'%s' must be specified either as a variable or in 'attr(x, \"%s\")'", n, a))
    }
    if (missing(site))
      site <- checkattr("name", "site")
    st[2] <- site
    if (missing(lat))
      lat <- checkattr("latitude", "lat")
    if (missing(lon))
      lon <- checkattr("longitude", "lon")
    if (missing(alt))
      alt <- checkattr("elevation", "alt")
    lla <- sprintf("\t%s", paste(lat, lon, alt, sep="\t"))  # leading tab
    st <- c("[SITE]", site,
            "[LAT, LON and ALT]", lla,
            "[WEATHER FILES]", basename(datfile),
            "[FORMAT]", NA, "[END]")
  }
  # make dat/sr file first
  if (!"date" %in% names(x) || !inherits(x$date, c("Date", "POSIXct")))
    stop("There must be a 'x$date' column with class 'Date' or 'POSIXct'")
  rng <- range(x$date, na.rm=TRUE)
  sdat <- data.frame(date=seq(rng[1], rng[2], by="day"))
  sdat$year <- as.integer(format(sdat$date, "%Y"))
  sdat$jday <- as.integer(format(sdat$date, "%j"))
  x$year <- x$yday <- x$jday <- NULL
  m <- match(names(x), conv$seas)
  m <- sort(m[!is.na(m)])
  x <- x[,c("date", conv$seas[m])]
  names(x) <- c("date", conv$lars[m])
  x <- merge(sdat, x, by="date", all.x=TRUE)
  x$date <- NULL
  names(x) <- toupper(names(x))
  write.table(x, datfile, quote=FALSE, row.names=FALSE, col.names=FALSE,
              na="-99", dec=".", sep="\t")
  if (mkst) {
    st[8] <- paste(names(x), collapse=" ")  # format
    writeLines(st, stfile)
  }
}

"write.lars.scenario" <-
function(file, x1, x2, name = "anomaly") {
  t_max <- change(x1, x2, "t_max")
  t_min <- change(x1, x2, "t_min")
  if (!"t_mean" %in% names(x1))
    x1$t_mean <- apply(x1[,c("t_min", "t_max")], 1, mean, na.rm=TRUE)
  if (!"t_mean" %in% names(x2))
    x2$t_mean <- apply(x2[,c("t_min", "t_max")], 1, mean, na.rm=TRUE)
  t_mean <- change(x1, x2, "t_mean")
  precip <- change(x1, x2, "precip", disc=TRUE, inter=TRUE)
  solar <- change(x1, x2, "solar")
  tb <- data.frame(row.names=month.abb, rain=rep(1, 12),
                   wet=1, dry=1, min=0, max=0, tsd=1, rad=1)
  tb$rain <- precip$cent.rel[,"precip"]
  tb$wet <- precip$wet.rel
  tb$dry <- precip$dry.rel
  tb$min <- t_min$cent.abs
  tb$max <- t_max$cent.abs
  tb$tsd <- t_mean$sprd.rel
  tb$rad <- solar$cent.rel
  if (file == "") {
    con <- stdout()
  } else {
    con <- file(file, "w")
    on.exit(close(con))
  }
  writeLines(c("// mon\train\twet\tdry\tmin\tmax\ttsd\trad",
               "[NAME]", name, "[GCM PREDICTIONS]"), con)
  write.table(round(tb, 4), con, quote=FALSE, sep="\t", col.names=FALSE)
  writeLines("[END]", con)
  invisible(tb)
}

"lars2help" <-
function(infile, outfile, year.offset, site) {
  fp <- file(infile, "r")
  l <- readLines(fp)
  close(fp)

  name <- l[2]
  lat <- as.numeric(substr(l[4], 1, 9))
  infile2 <- l[6]
  message("name: ", name, "\nlat: ", lat)

  message("reading ", infile2)
  x <- read.lars(infile2, year.offset=year.offset)
  x$t_mean <- rowMeans(x[,c("t_min", "t_max")])

  outfiles <- paste(outfile, c(4, 7, 13), sep=".D")
  names(outfiles) <- c("precip", "temp", "solar")

  message("writing ", outfiles["precip"])
  write.help(outfiles["precip"], x, "precip", name, site, lat)
  message("writing ", outfiles["temp"])
  write.help(outfiles["temp"], x, "t_mean", name, site, lat)
  message("writing ", outfiles["solar"])
  write.help(outfiles["solar"], x, "solar", name, site, lat)
  message("done")
}
