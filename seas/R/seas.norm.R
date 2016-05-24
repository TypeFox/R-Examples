"seas.norm" <-
function(x, var, fun = "median", norm = "days", year.filter,
         ann.only = FALSE, precip.norm = FALSE) {
  orig <- as.character(substitute(x))[[1]]
  if (precip.norm && (any(!c("rain", "snow") %in% x$var))) {
    warning(
      paste(
        gettextf("%s does not have either %s or %s variables, or both",
                 sQuote(orig), sQuote("rain"), sQuote("snow")),
            gettextf("using only %s", var), sep="\n"))
    precip.norm <- FALSE
  }
  x <- seas.sum.check(x, orig, var, norm, year.filter, ann.only)
  var <- x$var
  n.years <- length(x$years)
  if (is.function(fun))
    fun <- as.character(substitute(fun))[[1]]
  num <- length(x$bins)
  ann <- data.frame(var=NA)  # determine annual normals
  var.a <- x$ann[,var]
  if (any(is.finite(var.a))) {  # are values to calculate annual normals?
    ann$var <- eval(call(fun, var.a, na.rm=TRUE))
    if (precip.norm) {
      rain.a <- x$ann$rain
      snow.a <- x$ann$snow
      ann$rain <- eval(call(fun, rain.a, na.rm=TRUE))
      ann$snow <- eval(call(fun, snow.a, na.rm=TRUE))
    }
    ann$days <- eval(call(fun, x$ann$days, na.rm=TRUE))
    ann$active <- NA
    ann$na <- eval(call(fun, x$ann$na, na.rm=TRUE))
  } else {  # incomplete years
    ann$var <- NA
    if (precip.norm)
      ann$snow <- ann$rain <- rain.a <- snow.a <- NA
    ann$active <- ann$days <- NA
    ann$na <- NA
    warning(gettext("not enough data to determine annual normals"))
    if (fun == "median") {
      warning(gettextf("changing %s from %s to %s",
                       sQuote("fun"), sQuote(fun), sQuote("mean")))
      fun <- "mean"
    }
  }
  if (x$a.cut) {
    active.a <- x$ann$active
    ann$active <- if (any(is.finite(active.a)))
      eval(call(fun, active.a, na.rm=TRUE)) else 0
  } else {
    ann$active <- NULL
  }
  if (!ann.only) {
    ## .b suffix is a matrix of sums of bin for each year
    var.b <- x$seas[,,var] # drop=TRUE is default
    norm.b <- x$norm[,,var]
    var.n <- x$seas/x$norm # normalized
    if (dim(var.n)[1] == 1)
      var.n <- var.n[1,,]
    var.n[x$norm <= 0] <- 0
    if (precip.norm) {
      rain.b <- x$seas[,,"rain"]
      snow.b <- x$seas[,,"snow"]
    }
    seas <- data.frame(var=1:num * NA, row.names=x$bins)
    if (fun == "median" && n.years <= 2) {
      if (n.years != 1)
        warning("not enough years to find the median; using mean")
      fun <- "mean"
    }
    if (n.years == 1) {
      seas$var <- var.n[,var]
      if (precip.norm) {
        seas$rain <- var.n[, "rain"]
        seas$snow <- var.n[, "snow"]
      }
      if (x$a.cut)
        seas$active <- x$active[1,,var] / norm.b
      seas$days <- x$days[1,]
      seas$na <- x$na[1,] / norm.b
    } else if (fun == "median" && n.years > 2) {
      quan <- data.frame(var=NA)
      secant <- function (f) {
        ## Secant method to find a root; f is a function which needs to
        ## be zero, specifically a quantile in [0, 1]
        ## Cheny, E. W. and Kincaid, D. 1999
        a <- 0.45; b <- 0.55 # start around the 50% quantile or true median
        fa <- f(a); fb <- f(b)
        while (fa == fb && b < 1) # if there is no slope around initial guess
          fb <- f(b <- b + 0.05)
        for (i in 1:50) {
          if (a == b) stop("can not converge; outside quantile bounds [0, 1]")
          if (abs(fa) > abs(fb)) {  # swap
            z <- a;   a <- b;   b <- z
            fz <- fa; fa <- fb; fb <- fz
          }
          d <- (b - a) / (fb - fa)
          b <- a; fb <- fa
          d <- d * fa
          if (abs(fa) < 0.0001) return(a)  # normal exit point for function
          a <- a - d
          if (a < 0) a <- 0 # fix the bounds between [0, 1]
          if (a > 1) a <- 1
          fa <- f(a)
        }
        warning(gettextf("no convergence of quantile; %f > 0.0001", abs(fa)))
        return(a)
      }
      quan$var <- secant(function(qu)
                           return(ann$var - sum(apply(var.b, 2,
                                                      quantile, qu,
                                                      na.rm=TRUE,
                                                      names=FALSE))))
      seas$var <- apply(var.n[,,var], 2,
                          quantile, quan$var,
                          na.rm=TRUE, names=FALSE)
      if (precip.norm) {
        quan$snow <- quan$rain <- NA
        ## calculate the quantile of annual rain + snow, which equals
        ## the anual precipitation
        quan$rainsnow <- secant(function(qu)
                                return(ann$var
                                       - quantile(rain.a, qu,
                                                  na.rm=TRUE,
                                                  names=FALSE)
                                       - quantile(snow.a, qu,
                                                  na.rm=TRUE,
                                                  names=FALSE)))
        rain.am <- quantile(ann$rain, quan$rainsnow,
                            na.rm=TRUE, names=FALSE)
        ## calculate the quantile of all rain data in the bins to
        ## equal the annual rain volume
        quan$rain <- secant(function(qu)
                            return(rain.am - sum(apply(rain.b, 2,
                                                       quantile, qu,
                                                       na.rm=TRUE,
                                                       names=FALSE))))
        seas$rain <- apply(var.n[,,"rain"], 2,
                           quantile, quan$rain,
                           na.rm=TRUE, names=FALSE)
        snow.am <- quantile(ann$snow, quan$rainsnow,
                            na.rm=TRUE, names=FALSE)
        if (snow.am > 0) {  # median may be zero
          quan$snow <- secant(function(qu)
                              return(snow.am - sum(apply(snow.b, 2,
                                                         quantile, qu,
                                                         na.rm=TRUE,
                                                         names=FALSE))))
          seas$snow <- apply(var.n[,,"snow"], 2, quantile,
                             quan$snow,
                             na.rm=TRUE, names=FALSE)
        } else seas$snow <- seas$rain * 0
        ## calculate the fraction of rain
        rs.f <- with(seas, rain / (rain + snow))
        if (any(is.nan(rs.f)))
          rs.f[seas$rain + seas$snow == 0] <- 1  # avoid divide by zero
        seas$rain <- seas$var * rs.f
        seas$snow <- seas$var * (1 - rs.f)
        ## these two operations make the bars equal for
        ## precip.norm=TRUE and FALSE
      }
      if (x$a.cut && ann$active > 0) {
        quan$active <- secant(function(qu)
                              return(ann$active - sum(apply(x$active, 2,
                                                            quantile, qu,
                                                            na.rm=TRUE,
                                                            names=FALSE))))
        seas$active <- apply(x$active[,,var] / norm.b, 2,
                             quantile, quan$active,
                             na.rm=TRUE, names=FALSE)
      }
      if (all(apply(x$days, 1, function(r)(r == x$days[1,]))))
        quan$days <- 0.5  # secant method fails if all same value
      else
        quan$days <- secant(function(qu)
                            return(ann$days - sum(apply(x$days, 2,
                                                        quantile, qu,
                                                        na.rm=TRUE,
                                                        names=FALSE))))
      seas$days <- apply(x$days, 2, quantile, quan$days,
                         na.rm=TRUE, names=FALSE)
      if (ann$na > 0) {
        quan$na <- secant(function(qu)
                          return(ann$na - sum(apply(x$na, 2,
                                                    quantile, qu,
                                                    na.rm=TRUE,
                                                    names=FALSE))))
        seas$na <- apply(x$na / norm.b, 2, quantile, quan$active,
                         na.rm=TRUE, names=FALSE)
      } else seas$na <- quan$days * 0
    } else {  # calculate stats conventionally (much faster)
      seas$var <- apply(var.n[,,var], 2, fun, na.rm=TRUE)
      if (precip.norm) {
         seas$rain <- apply(var.n[,,"rain"], 2, fun, na.rm=TRUE)
         seas$snow <- apply(var.n[,,"snow"], 2, fun, na.rm=TRUE)
      }
      if (x$a.cut)
        seas$active <- apply(x$active[,,var] / norm.b, 2, fun, na.rm=TRUE)
      seas$days <- apply(x$days, 2, fun, na.rm=TRUE)
      seas$na <- apply(x$na / norm.b, 2, mean, na.rm=TRUE)
    }
  }  # end if (!ann.only)
  fixnames <- function(n) {
    n[n %in% "var"] <- var
    return(n)
  }
  names(ann) <- fixnames(names(ann))
  l <- list(ann=ann)
  if (!ann.only) {
    attr(l, "class") <- "seas.norm"
    names(seas) <- fixnames(names(seas))
    l$seas <- seas
    l$width <- x$width
    l$bins <- x$bins
    l$bin.lengths <- x$bin.lengths
  }
  l$start.day <- x$start.day
  l$year.range <- x$year.range
  l$var <- var
  l$long.name <- x$long.name[[var]]
  l$units <- x$units[[var]]
  l$norm <- norm
  l$ann.only <- ann.only
  l$precip.norm <- precip.norm
  l$a.cut <- x$a.cut
  if (n.years > 1)
    l$fun <- fun
  if (!ann.only && fun == "median") {
    names(quan) <- fixnames(names(quan))
    l$quantile <- quan
  }
  l$id <- x$id
  l$name <- x$name
  l
}

"precip.norm" <-
function(x, fun = "median", norm = "days", year.filter) {
  if (is.function(fun))
    fun <- as.character(substitute(fun))[[1]]
  seas.norm(x, norm=norm, year.filter=year.filter, fun=fun, precip.norm=TRUE)
}
