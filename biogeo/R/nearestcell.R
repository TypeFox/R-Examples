nearestcell <-
function (dat, rst) 
{
  length(dat)
  length(rst)
  fieldsmissing(dat, fields = c("ID", "x", "y", "x_original", 
                                "y_original", "Correction", "Modified", "Exclude","Reason"))
  dat$Correction <- as.character(dat$Correction)
  fx <- which(dat$Exclude == 0)
  x1 <- coord2numeric(dat$x[fx])
  y1 <- coord2numeric(dat$y[fx])
  datid <- dat$ID[fx]
  dd <- data.frame(x1, y1)
  ce0 <- cellFromXY(rst, dd)
  vals <- extract(rst, dd)
  f <- which(is.na(vals))
  if(length(f)==0)
    stop('There are no missing values')
  id <- datid[f]
  ce1 <- cellFromXY(rst, dd[f, ])
  ff <- which(!is.na(ce1))
  ce3 <- ce1[ff]
  dd2 <- dd[ff, ]
  id2 <- id[ff]
  bb <- {
  }
  for (i in 1:length(ce3)) {
    a <- adjacent(rst, ce3[i], directions = 8, pairs = FALSE, 
                  target = NULL, sorted = FALSE, include = FALSE, id = FALSE)
    idx <- id2[i]
    b <- data.frame(i, a, id2 = idx)
    bb <- rbind(bb, b)
  }
  xy <- xyFromCell(rst, bb$a, spatial = FALSE)
  vals <- extract(rst, xy)
  g <- data.frame(bb, xy, vals)
  if (any(!is.na(vals))) {
    fv <- which(!is.na(vals))
    id <- bb$i[fv]
    uid <- unique(id)
    g1 <- na.omit(g)
    near <- {
    }
    for (j in 1:length(uid)) {
      uj <- uid[j]
      fx <- which(g1$i == uj)
      nce <- g1[fx, ]
      nce2 <- as.matrix(cbind(nce$x, nce$y))
      pts <- dd2[uj, ]
      dst <- pointDistance(c(pts$x1, pts$y1), nce2, longlat = F)
      fm <- which.min(dst)
      nr <- nce[fm, ]
      near <- rbind(near, nr)
    }
  } else{
    stop('There are no records close enough to the nearest land/sea cells')
  }
  mod <- format(Sys.time(), "%d-%m-%Y %H:%M:%S")
  dat$Modified <- as.character(dat$Modified)
  for (i in 1:nrow(near)) {
    f <- which(dat$ID == near$id2[i])
    dat$x_original[f] <- dat$x[f]
    dat$y_original[f] <- dat$y[f]
    dat$x[f] <- near$x[i]
    dat$y[f] <- near$y[i]
    dc <- dat$Correction[f]
    dc <- str_replace(dc, "[.]", "7")
    dat$Correction[f] <- dc
    dat$Modified[f] <- mod
  }
  moved <- data.frame(ID = near$id2, x = near$x, y = near$y)
  datx <- list(dat = dat, moved = moved)
  return(datx)
}
