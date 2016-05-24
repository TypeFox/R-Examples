keypoints <- function(x, ...) {
  UseMethod("keypoints")
}

# tracks

keypoints.trkset <- function(x, ...) {
  val <- NULL
  for(i in 1:length(x))
    val <- rbind(val, keypoints.trk(x[[i]]))
  return(val)
}


keypoints.trk <- function(x, ...) {
  # row of the track
  nrow <- dim(x$data)[1]
  
  # Te shift for BTO
  shiftTe <- 100
  
  zams <- 171
  exhc <- 571
  flash <- 871

  # search TO
  nTO <- which(x$data$logTe == max(x$data$logTe))
  TO <- sapply(x$data[nTO,], mean)

  # search BTO
  sel <- which(x$data$logTe < log10(10^(TO["logTe"]) - shiftTe))
  nBTO <- sel[sel > nTO[1]][1]
  if(!is.na(nBTO))  
    BTO <- x$data[nBTO,]

  # build the output data frame 
  val <- x$data[c(zams),]
  val <- rbind(val, TO)

  # optional BTO point (it does not exist for all tracks)
  if(!is.na(nBTO))  
    val <- rbind(val, BTO)

  val <- rbind(val, x$data[exhc,])

  # points for late evolutionary phases (not for all tracks)
  if(nrow > exhc) {
    val <- rbind(val, x$data[flash,])
    rownames(val) <- c("ZAMS", "TO", "BTO", "exHc", "Heflash")
  } else {
    if(!is.na(nBTO))  
      rownames(val) <- c("ZAMS", "TO", "BTO", "exHc")
    else
      rownames(val) <- c("ZAMS", "TO", "exHc")
  }

                                        # set Hc = 0
                                        # (otherwise it shows the He abund)
  val["exHc","Hc"] <- 0
  if(nTO[1] == exhc)
     val["TO","Hc"] <- 0
  
                                        # add track characteristics
  val$M <- x$mass
  val$z <- x$z
  val$y <- x$y
  val$ml <- x$ml
  val$alpha.enh <- x$alpha.enh

                                        # add phase identifier
                                        # 1=ZAMS, 2=TO, 3=BTO, 4=exHc,
                                        # 5=He flash
  if(is.na(nBTO))
    val$id <- c(1,2,4)
  else
    val$id <- seq(1,dim(val)[1], by=1)
    
  return(val)
}

# isochrones

keypoints.isoset <- function(x, ...) {
  val <- NULL
  for(i in 1:length(x))
    val <- rbind(val, keypoints.iso(x[[i]]))
  return(val)
}

keypoints.iso <- function(x, ...) {

                                        # Te shift for BTO
  shiftTe <- 100

                                        # search TO
  nTO <- which(x$data$logTe == max(x$data$logTe))
  TO <- sapply(x$data[nTO,], mean)
  
                                        # search BTO
  nBTO <- which(x$data$logTe < log10(10^(TO["logTe"]) - shiftTe) & x$data$logL > TO["logL"])[1]
  if(!is.na(nBTO))  
    BTO <- x$data[nBTO,]

  val <- NULL
  val <- rbind(val, TO)
  val <- rbind(val, BTO)
  
                                        # add iso characteristics
  val$age <- x$age
  val$z <- x$z
  val$y <- x$y
  val$ml <- x$ml
  val$alpha.enh <- x$alpha.enh

  
                                        # add phase identifier
                                        # 1=TO, 2=BTO
  if(!is.na(nBTO)) {
    val$id <- c(1,2)
    rownames(val) <- c("TO", "BTO")
  } else {
    val$id <- 1
    rownames(val) <- c("TO")
  }
  
  return(val)
}
