## General functions useful for processing the Axion data.
## 2013-01-29

## This variable stores all the information related to the wells; typically this
## is accessed through the plateinfo(arrayname) function.

.plateinfo <- list("Axion 48 well"=list(
                     n.well=48,
                     wells=paste( rep(LETTERS[6:1], each=8), rep(1:8,6), sep=''),
                     n.well.r=6,
                     n.well.c=8,
                     layout=c(8,6),
                     n.elec.r=4,
                     n.elec.c=4),
                 "Axion 12 well"=list(
                     n.well=12,
                     wells=paste( rep(LETTERS[3:1], each=4), rep(1:4,3), sep=''),
                     n.well.r=3,
                     n.well.c=4,
                     layout=c(4,3),
                     n.elec.r=8,
                     n.elec.c=8))

                     
plateinfo <- function(arrayname) {
    ## Return useful information related to arrayname
    ## 
    ## plateinfo("Axion 12 well")
    res <- .plateinfo[[arrayname]]
    if (is.null(res)) {
        stop("arrayname not recognised:", arrayname)
    } else {
        res
    }
}
    
## * Scripts to convert data into HDF5.
axion.elec.name.to.xy <- function(name, plateinfo) {
  ## Convert electrode name to  (x,y) position.
  ## plateinfo stores all the information about the plates.
  ## and hence the well layout of the plate.

  max.well.row <-  plateinfo$n.well.r
  max.well.col <-  plateinfo$n.well.c
  max.elec.row <-  plateinfo$n.elec.r
  max.elec.col <-  plateinfo$n.elec.c
  
  
  well.r <- max.well.row - match(substring(name, 1,1), LETTERS)
  well.c <- as.integer(substring(name, 2,2)) - 1
  elec.r <- as.integer(substring(name, 5,5)) - 1
  elec.c <- as.integer(substring(name, 4,4)) - 1

  gap <- 1
  spacing <- 200                        #electrode spacing.
  well.wid <- (max.elec.col+gap)*spacing
  well.ht  <- (max.elec.row+gap)*spacing

  x <- (well.c*well.wid) + (elec.c*spacing)
  y <- (well.r*well.ht)  + (elec.r*spacing)

  cbind(x,y)

}

axion.spikesum <- function(spikes) {
  ## Generate a simple summary of the spikes list.
  len <- length(spikes)
  all.range <- sapply(spikes, range)
  nspikes <- sum(sapply(spikes, length))
  min <- min(all.range[1,])
  max <- max(all.range[2,])
  str <- sprintf("summary: %d electrodes %d spikes, min %.4f max %.4f",
                 len, nspikes, min, max)
  str
}

axion.spikesum2 <- function(spikes) {
  ## Generate a simple summary of the spikes list.
  ## This version returns a vector, rather than a string.  This is more
  ## useful for building a data frame of values.
  len <- length(spikes)
  all.range <- sapply(spikes, range)
  nspikes <- sum(sapply(spikes, length))
  min <- min(all.range[1,])
  max <- max(all.range[2,])
  str <- sprintf("summary: %d electrodes %d spikes, min %.4f max %.4f",
                 len, nspikes, min, max)
  ##str
  c(nelectrodes=len, nspikes=nspikes, time.min=min, time.max=max)
}

axion.spikestodf <- function(spikes) {
  ## Convert a list of spikes to a 2-column  (elec, time) data frame.
  names <- names(spikes)
  names(spikes) <- NULL
  nspikes <- sapply(spikes, length)
  data.frame(elec=rep(names, times=nspikes), time=unlist(spikes))
}


axion.guess.well.number <- function(channels) {
  ## Given the channel names, guess the number of wells on the plate.
  ## This works on the logic that certain electrode names will only be
  ## found on certain plates. e.g. the electrode name "D6_33" can only appear
  ## on a well with 48 arrays.
  ##
  ## axion.guess.well.number("D3_33")  ## should be 48.
  ## axion.guess.well.number("B3_53")  ## should be 12
  ## axion.guess.well.number("A2_11") ## this is ambiguous.
  
  well.r <- match(substring(channels, 1,1), LETTERS)
  well.c <- as.integer(substring(channels, 2,2))
  elec.r <- as.integer(substring(channels, 5,5))
  elec.c <- as.integer(substring(channels, 4,4))

  max.well.r <- max(well.r)
  max.well.c <- max(well.c)

  max.elec.r <- max(elec.r)
  max.elec.c <- max(elec.c)

  nplates <- length(.plateinfo)
  well <- 0
  for (i in 1:nplates) {
    plateinfo <- .plateinfo[[i]]
    if (max.well.r <= plateinfo$n.well.r &&
        max.well.c <= plateinfo$n.well.c &&
        max.elec.r <= plateinfo$n.elec.r &&
        max.elec.c <= plateinfo$n.elec.c) {
      well <- plateinfo$n.well
      break;
    }
  }
  if (well == 0) {
    stop("Cannot guess number of wells on plate.")
  }
  well
}


axion.electrodes.on.well <- function(well, electrodes) {
  ## Return names of electrodes that are on well WELL.
  matches <- grep(well, electrodes)
  electrodes[matches]
}

axion.elec2well <- function(elec) {
  ## Extract well name from ELECtrode name.
  substring(elec, 1, 2)
}


#this function adds well info and nAE=# active electrodes
.well.info<-function(s2) {
  if ("12"==substring(s2$layout$array,7,8)){
    s2$wells<-c("A1","A2","A3","A4","B1","B2","B3","B4","C1","C2","C3","C4")
  } else if ("48"==substring(s2$layout$array,7,8)){
    s2$wells<-c('A1','A2','A3','A4','A5','A6','A7','A8',
                'B1','B2','B3','B4','B5','B6','B7','B8',
                'C1','C2','C3','C4','C5','C6','C7','C8',
                'D1','D2','D3','D4','D5','D6','D7','D8',
                'E1','E2','E3','E4','E5','E6','E7','E8',
                'F1','F2','F3','F4','F5','F6','F7','F8')
  }
  #add number of active electrodes
  s2$nAE<-rep(0,length(s2$wells))
  names(s2$nAE)<-s2$wells
  for (i in 1:s2$NCells){
    s2$nAE[which(substr(s2$channels[i],1,2)==(s2$wells))]=
      s2$nAE[which(substr(s2$channels[i],1,2)==(s2$wells))]+1
    s2$cw[i]<-substr(s2$channels[i],1,2)
  }
  s2
}
