faremalm2 <- function(dat = NULL, noutput = 1, id = "id", year =
                      "year"){
  require(gdata)
  require(Hmisc)
  require(psych)

  dat <- dat[order(dat[[id]], dat[[year]]),]
  ids <- as.character(unique(dat[[id]]))
  years <- sort(unique(dat[[year]]))

  xy.idx <- 3:ncol(dat)

  Dt2t2 <- list()
  Dtt2 <- list()
  Dt2t <- list()

  cat("Dt2t2\n")
  for(j in years){
    tmp <- dat[dat[[year]] ==j,]
    base <- tmp[, xy.idx]
    frontier <- tmp[, xy.idx]
    if(nrow(tmp) > 100)
      Dt2t2[[paste("X",j, sep = "")]] <-
        unlist(1/effdea.b.f(base, frontier, noutput = noutput,
                            orientation = 2, rts = 1))
    else
      Dt2t2[[paste("X",j, sep = "")]] <-
        1/dea(base, frontier, noutput = noutput, orientation = 2,
                  rts = 1)[,1]
    names(Dt2t2[[paste("X",j, sep = "")]]) <-
      paste("id", tmp[[id]], sep = "")
  }

  Dt2t2 <- data.frame(Dt2t2 = unlist(Dt2t2))
  row.names(dat) <- paste("X",  dat[[year]], ".id",
                          dat[[id]], sep = "")
  dat <- merge(dat, Dt2t2, by.x = "row.names", by.y = "row.names",
  all.x = TRUE)
  row.names(dat) <- dat$Row.names
  dat$Row.names <- NULL

  # Dtt2
  cat("Dtt2\n")
  for(j in years[-1]){
    tmp.b <- dat[dat[[year]] ==j,]
    tmp.f <- dat[dat[[year]] ==j-1,]
    base <- tmp.b[, xy.idx]
    frontier <- tmp.f[, xy.idx]
    if(nrow(tmp.b) > 100)
      Dtt2[[paste("X",j, sep = "")]] <-
        unlist(1/effdea.b.f(base, frontier, noutput = noutput,
                            orientation = 2, rts = 1))
    else
      Dtt2[[paste("X",j, sep = "")]] <-
        1/dea(base, frontier, noutput = noutput, orientation = 2,
                  rts = 1)[,1]
    names(Dtt2[[paste("X",j, sep = "")]]) <-
      paste("id", tmp.b[[id]], sep = "")
  }

  Dtt2 <- data.frame(Dtt2 = unlist(Dtt2))
  row.names(dat) <- paste("X",  dat[[year]], ".id",
                          dat[[id]], sep = "")
  dat <- merge(dat, Dtt2, by.x = "row.names", by.y = "row.names",
  all.x = TRUE)
  row.names(dat) <- dat$Row.names
  dat$Row.names <- NULL


  # Dt2t
  cat("Dt2t\n")
  for(j in years[-length(years)]){
    tmp.f <- dat[dat[[year]] ==j+1,]
    tmp.b <- dat[dat[[year]] ==j,]
    base <- tmp.b[, xy.idx]
    frontier <- tmp.f[, xy.idx]
    if(nrow(tmp.f) > 100)
      Dt2t[[paste("X",j+1, sep = "")]] <-
        unlist(1/effdea.b.f(base, frontier, noutput = noutput,
                            orientation = 2, rts = 1))
    else
      Dt2t[[paste("X",j+1, sep = "")]] <-
        1/dea(base, frontier, noutput = noutput, orientation = 2,
                  rts = 1)[,1]
    names(Dt2t[[paste("X",j+1, sep = "")]]) <-
      paste("id", tmp.b[[id]], sep = "")
  }

  Dt2t <- data.frame(Dt2t = unlist(Dt2t))
  row.names(dat) <- paste("X",  dat[[year]], ".id",
                          dat[[id]], sep = "")
  dat <- merge(dat, Dt2t, by.x = "row.names", by.y = "row.names",
               all.x= TRUE)
  row.names(dat) <- dat$Row.names
  dat$Row.names <- NULL

  dat <- dat[order(dat[[id]], dat[[year]]),]
  dat$Dtt <- unlist(by(dat$Dt2t2, dat[[id]], Lag))

  dat <- dat[order(dat[[year]], dat[[id]]),]
  dat <- na.omit(dat)

  dat$ec <- dat$Dt2t2/dat$Dtt
  dat$tc <- ((dat$Dtt2/dat$Dt2t2)*(dat$Dtt/dat$Dt2t))^.5
  dat$pc <- dat$ec * dat$tc
  
  return(dat)

}
