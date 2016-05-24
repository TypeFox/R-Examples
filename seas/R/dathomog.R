"dathomog" <-
function (x1, x2, by="date", plot=FALSE) {
  orig1 <- as.character(substitute(x1))[[1]]
  orig2 <- as.character(substitute(x2))[[1]]
  sc1 <- seas.df.check(x1, orig1)
  sc2 <- seas.df.check(x2, orig2)
  n1 <- names(x1)
  n2 <- names(x2)
  a <- list() # copy attributes
  for (n in n2)
    a[[n]] <- attributes(x2[[n]])
  for (n in n1)
    a[[n]] <- attributes(x1[[n]])
  if (!(by %in% n1 && by %in% n2))
    stop(gettextf("problems encountered while trying to find %s in names of 'x1' and 'x2' for %s argument\n",
                  sQuote(by),sQuote("by")))
  vars.all <- union(n1, n2)  # keep `by' in this one
  vars <- vars.all[!vars.all %in% by]  # stip `by' out
  vars.x <- paste(vars[vars %in% n1], ".x", sep="")
  vars.y <- paste(vars[vars %in% n2], ".y", sep="")
  dat <- merge(x1, x2, by=by, all=TRUE)
  if (plot) {
    require(MASS)
    par(ask=TRUE)
    var <- 1:length(vars)
    var <- var[sapply(dat[,vars.x], is.numeric)]
    RsqrSym <- paste("R", iconv("\262", "latin1", ""), sep="")
    for (v in var) {
      var <- vars[v]
      x <- as.numeric(dat[,vars.x[v]])
      y <- as.numeric(dat[,vars.y[v]])
      cr <- !is.na(x) & !is.na(y)
      x <- x[cr]; y <- y[cr]
      n <- length(x)
      suppressWarnings(lim <- range(x, y))
      if (n < 1 || diff(lim) == 0) {
        frame()
        title(var, xlab=n1, ylab=n2)
        if (n < 1)
          text(0.5, 0.5, gettext("no overlap of data"))
        else
          text(0.5, 0.5, paste(gettext("no variance to data"),
                           sprintf("n = %i; value = %f", x[1]), sep="\n"))
      } else {  # calculate least-squares perpendicular offsets line
        B <- 0.5 * ((sum(y^2) - n * mean(y)^2) - (sum(x^2) - n * mean(x)^2))/
          (n * mean(x) * mean(y) - sum(x * y))
        slope <- (-B + sqrt(B^2 + 1))
        inter <- mean(y) - slope*mean(x)
        rsq <- cor(x, y)^2
        lim <- lim + diff(lim) * (4 * c(-1, 1) / 100)  # expand range +/-4%
        image(kde2d(x, y, c(3, 3), 50, c(lim, lim)),
              col=rev(terrain.colors(30)), xlim=lim, ylim=lim,
              xlab=n1, ylab=n2, main=var)
        if (n > 200) { x <- jitter(x); y <- jitter(y) }
        if (n < 1000) points(x, y)
        abline(inter, slope)
        legend(lim[1], lim[2],
               paste(c("inter. ", "slope", paste(RsqrSym, "  "), "n       "),
                     round(c(inter, slope, rsq, n), 3), sep=": "), bg="white")
      }
    }
  }

  by.c <- dat[,by]
  dat.x <- dat[,vars.x]
  dat.y <- dat[,vars.y]
  dat <- data.frame(by=by.c)
  dat$by <- NULL
  dat[,vars.all] <- NA
  dat[,by] <- by.c
  names(dat.x) <- vars
  dat[,vars] <- dat.x
  for (v in names(dat)) {
    dat.na <- is.na(dat[[v]])
    dat[[v]][dat.na] <- dat.y[[paste(v, ".y", sep="")]][dat.na]
  }
  for (n in names(dat))
    attributes(dat[[n]]) <- a[[n]]
  invisible(dat)
}
