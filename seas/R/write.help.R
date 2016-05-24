"write.help" <-
function(file, dat, var, name = "", region = "", lat,
         visual.help = FALSE, metric = TRUE) {
  if (!var %in% c("t_mean", "precip", "solar"))
    stop(gettextf("%s must be one of %s, %s or %s for daily mean temperature, precipitation and global solar radiation",
                  sQuote("var"), sQuote("t_mean"),
                  sQuote("precip"), sQuote("solar")))
  orig <- as.character(substitute(dat))[[1]]
  sc <- seas.df.check(dat, orig, var)
  ## The header for HELP is as follows:
  ## Line 1: %2i # (1) data source flag; use 3 for "create/edit"
  ## Line 2: %2i # (2) units; imperial=1, metric=2
  ## Line 3: %-20s %-20s # (3) name and region
  ## Line 4: depends on file type
  header <- formatC(c(3, ifelse(metric, 2, 1), NA, NA), 0, 2)
  header[3] <- sprintf("%-20s %-20s", name, region)
  if (var == "solar") {  # solar files require latitude
    header[4] <- if (visual.help)
      sprintf("%6.2f", lat)
    else
      sprintf("%9.2f", lat)
    type <- 5
  } else if (var %in% c("precip", "t_mean")) {  # these file require normals
    if (var == "precip") {
      seas.nm <- seas.norm(seas.sum(dat, var="precip", width="mon"),
                           fun=mean)$seas
      nm <- seas.nm$precip * seas.nm$days  # in mm/month
      type <- 1
    } else {  # t_mean
      nm <- tapply(dat$t_mean, mkseas(dat, width="mon"), mean, na.rm=TRUE)
      type <- 3
    }
    header[4] <- paste(sprintf("%6.1f", nm), collapse="")
  }
  header <- paste(header, collapse="\n")
  if (visual.help)
    type <- type + 1
  yearrange <- range(as.integer(format(dat$date, "%Y")))
  years <- seq(yearrange[1], yearrange[2], by=1)
  days <- rep(365, length(years))
  days[years %% 4 == 0 & years %% 100 != 0 | years %% 400 == 0] <- 366
  val <- dat[[var]]
  if (sum(days) != length(val))
    stop("Number of values should correspond to the number of possible days between whole years")
  if (any(!is.finite(val)))
    stop("Missing values are not allowed")
  invisible(.C("writeHELP", file, header, as.integer(type),
               as.integer(yearrange[1]), as.integer(diff(yearrange) + 1),
               as.double(val), PACKAGE="seas"))
}
