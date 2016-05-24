dmsparse <-
function (dat, x = "long", y = "lat", id = "ID") 
{
  cn <- names(dat)
  f1 <- match(x, cn)
  f2 <- match(y, cn)
  f3 <- match(id, cn)
  if (is.na(f1)) 
    stop(paste("The field", x, "does not exist"))
  if (is.na(f2)) 
    stop(paste("The field", y, "does not exist"))
  if (is.na(f3)) 
    stop(paste("The field", id, "does not exist"))
  r <- c(cn[f1], cn[f2], cn[f3])
  cnrem <- setdiff(cn, r)
  frem <- match(cnrem, cn)
  ident <- dat[, f3]
  x1 <- dat[, f1]
  y1 <- dat[, f2]
  nr <- nrow(dat)
  datx <- data.frame(ident, x1, xdeg = NA, xmin = NA, xsec = NA, 
                     EW = NA, x = NA, xnotes = "...", exclude = 0, stringsAsFactors = F)
  daty <- data.frame(ident, y1, ydeg = NA, ymin = NA, ysec = NA, 
                     NS = NA, y = NA, ynotes = "...", exclude = 0, stringsAsFactors = F)
  ftx <- fmtcheck(datx$x1)
  fx1 <- which(ftx == 0)
  fx0 <- which(ftx > 0)
  datx$xnotes[fx0] <- "incorrect format"
  datx$exclude[fx0] <- 1
  fty <- fmtcheck(daty$y1)
  fy1 <- which(fty == 0)
  fy0 <- which(fty > 0)
  daty$ynotes[fy0] <- "incorrect format"
  daty$exclude[fy0] <- 1
  suppressWarnings(xa <- as.numeric(datx$x1))
  suppressWarnings(ya <- as.numeric(daty$y1))
  fdx <- which(!is.na(xa))
  fdy <- which(!is.na(ya))
  datx$x[fdx] <- xa[fdx]
  daty$y[fdy] <- ya[fdy]
  ftx[fdx] <- 1
  fty[fdy] <- 1
  fx2 <- which(ftx == 0)
  fy2 <- which(fty == 0)
  dLx <- getletter(datx$x1[fx2])
  nax <- which(is.na(dLx$L))
  dLx$L[nax] <- "E"
  dLy <- getletter(daty$y1[fy2])
  nay <- which(is.na(dLy$L))
  dLy$L[nay] <- "N"
  dmsx <- sep(dLx)
  dmsy <- sep(dLy)
  ddx <- dms2dd(dmsx$deg, dmsx$min, dmsx$sec, dmsx$L)
  ddy <- dms2dd(dmsy$deg, dmsy$min, dmsy$sec, dmsy$L)
  dmsx$L[nax] <- NA
  ddx[nax] <- NA
  dmsy$L[nay] <- NA
  ddy[nay] <- NA
  datx$xdeg[fx2] <- dmsx$deg
  datx$xmin[fx2] <- dmsx$min
  datx$xsec[fx2] <- dmsx$sec
  datx$EW[fx2] <- dmsx$L
  daty$ydeg[fy2] <- dmsy$deg
  daty$ymin[fy2] <- dmsy$min
  daty$ysec[fy2] <- dmsy$sec
  daty$NS[fy2] <- dmsy$L
  datx$x[fx2] <- ddx
  daty$y[fy2] <- ddy
  x1 <- (abs(datx$x) > 180) * 1
  y1 <- (abs(daty$y) > 90) * 1
  if (any(na.omit(x1) == 1)) {
    fx1 <- which(x1 == 1)
    datx$x[fx1] <- NA
    datx$xnotes[fx1] <- "impossible coord value"
  }
  if (any(na.omit(y1) == 1)) {
    fy1 <- which(y1 == 1)
    daty$y[fy1] <- NA
    daty$ynotes[fy1] <- "impossible coord value"
  }
  dg1 <- (str_detect(datx$EW, "N")) * 1
  dg2 <- (str_detect(datx$EW, "S")) * 1
  dg <- na.omit(dg1 + dg2)
  if (any(dg >= 1)) {
    ff <- which((dg1 + dg2) >= 1)
    datx$xnotes[ff] <- "incorrect letter"
    datx$x[ff] <- NA
  }
  dg1 <- (str_detect(daty$NS, "E")) * 1
  dg2 <- (str_detect(daty$NS, "W")) * 1
  dg <- na.omit(dg1 + dg2)
  if (any(dg >= 1)) {
    ff <- which(c(dg1 + dg2) >= 1)
    daty$ynotes[ff] <- "incorrect letter"
    daty$y[ff] <- NA
  }
  fncx <- datx$x1%in%c("",NA)
  fncy <- daty$y1%in%c("",NA)
  datx$xnotes[fncx] <- "no coordinates"
  daty$ynotes[fncy] <- "no coordinates"
  exc <- (is.na(cbind(datx$x, daty$y))) * 1
  exc1 <- rowSums(exc)
  exclude <- ((datx$exclude + daty$exclude + exc1) >= 1) * 
    1
  datx <- datx[, -9]
  daty <- daty[, -c(1, 9)]
  data1 <- data.frame(datx, daty, dat[, frem], exclude, stringsAsFactors = F)
  names(data1) <- c("ID", "x_dms", "xdeg", "xmin", "xsec", 
                    "EW", "x", "xnotes", "y_dms", "ydeg", "ymin", "ysec", 
                    "NS", "y", "ynotes", cnrem, "exclude")
  return(data1)
}
