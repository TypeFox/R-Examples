#nts tiles, translated from Java

nts.MAP_SERIES_N_OF_80 <- matrix(c("910", "780", "560", "340", "120", 
                               NA, "781", "561", "341", "121"), ncol=5, byrow=TRUE)


nts.MAP_250K <- matrix(c("D", "C", "B", "A", 
                     "E", "F", "G", "H", 
                     "L", "K", "J", "I",
                     "M", "N", "O", "P"), ncol=4, byrow=TRUE)


nts.MAP_250K_N_OF_68 <- matrix(c("B", "A", "C", "D", "F", "E", "G", "H"), ncol=2, byrow=TRUE)


nts.MAP_50K <- matrix(c("04", "03", "02", "01", 
                    "05", "06", "07", "08",
                    "12", "11", "10", "09",
                    "13", "14", "15", "16"), ncol=4, byrow=TRUE)


nts.indexxy <- function(value, map) {
  which(map==value, arr.ind=TRUE)
}

nts.makebbox <- function(n, e, s, w) {
#   min max
#   x -66 -64
#   y  45  46
  matrix(c(w, s, e, n), byrow=FALSE, ncol=2, dimnames=list(c("x", "y"), c("min", "max")))
}


nts.widthandoffset250 <- function(lat) {
  if(lat >= 80.0) {
    c(8.0, 8.0)
  } else if(lat >= 68.0) {
    c(4.0, 0.0)
  } else {
    c(2.0, 0.0)
  }
}

nts.widthandoffsetseries <- function(lat) {
  if(lat >= 80) {
    c(16.0, 8.0)
  } else {
    c(8.0, 0.0)
  }
}


nts.mapsperseries <- function(tile250ky) {
  if(tile250ky >= 28) {
    2
  } else {
    4
  }
}


nts.tileseriesy <- function(lat) {
  as.integer(floor((lat-40.0) / 4.0))
}

nts.tileseriesx <- function(lon, lat) {
  wo <- nts.widthandoffsetseries(lat)
  as.integer(floor((lon+(144.0-wo[2]))/wo[1]))
}


nts.tileseries <- function(lon, lat) {
  c(nts.tileseriesx(lon, lat), nts.tileseriesy(lat))
}

nts.bboxseries <- function(tile) {
  minlat <- tile[2] * 4.0 + 40.0
  wo <- nts.widthandoffsetseries(minlat)
  minlon <- -144 + tile[1] * wo[1] + wo[2]
  nts.makebbox(minlat+4.0, minlon+wo[1], minlat, minlon) #n, e, s, w
}

nts.validtileseries <- function(tile) {
  if(tile[2] == 0) {
    if(7 <= tile[1] && tile[1] <= 11) {
      return(TRUE)
    }
  } else if(tile[2] == 1) {
    if(6 <= tile[1] && tile[1] <= 11 || tile[1] == 1 || tile[1] == 2) {
      return(TRUE)
    }
  } else if(2 <= tile[2] && tile[2] <= 4) {
    if(0 <= tile[1] && tile[1] <= 11) {
      return(TRUE)
    }
  } else if(5 <= tile[2] && tile[2] <= 8) {
    if(0 <= tile[1] && tile[1] <= 10) {
      return(TRUE)
    }
  } else if(tile[2] == 9) {
    if(0 <= tile[1] && tile[1] <= 9) {
      return(TRUE)
    }
  } else if(tile[2] == 10) {
    if(0 <= tile[1] && tile[1] <= 4) {
      return(TRUE)
    }
  } else if(tile[2] == 11) {
    if(1 <= tile[1] && tile[1] <= 4) {
      return(TRUE)
    }
  }
  FALSE
}

nts.idseries <- function(tile) {
  if(!nts.validtileseries(tile)) return(NA)
  
  if(tile[2]>=10) {
    return(nts.MAP_SERIES_N_OF_80[tile[2]-10+1, tile[1]+1])
  } else {
    seriesrow <- as.character(tile[2])
    seriescol <- as.character(11-tile[1])
    if(nchar(seriescol) == 1) {
      seriescol <- paste0("0", seriescol)
    }
    paste0(seriescol, seriesrow)
  }
}


nts.tileseriesbyid <- function(series) {
  if(nchar(series) >= 2) {
    result <- nts.indexxy(series, nts.MAP_SERIES_N_OF_80) #row, col returned
    if(length(result) > 0) {
      return(c(result[2]-1, result[1]+10-1))
    } else {
      seriesy <- as.integer(substr(series, nchar(series), nchar(series))) #last character
      seriesx <- 11- as.integer(substr(series, 1, nchar(series)-1))
      return(c(seriesx, seriesy))
    }
  } else {
    NULL
  }
}

nts.tileseriesfromtile250 <- function(tile250) {
  mapsperseries <- nts.mapsperseries(tile250[2])
  tilex <- as.integer(floor(tile250[1] / mapsperseries))
  tiley <- as.integer(floor(tile250[2] / 4))
  c(tilex, tiley)
}


nts.tile250y <- function(lat) {
  as.integer(floor(lat-40.0))
}

nts.tile250x <- function(lon, lat) {
  wo <- nts.widthandoffset250(lat)
  as.integer(floor((lon+(144-wo[2]))/wo[1]))
}

nts.tile250 <- function(lon, lat) {
  c(nts.tile250x(lon, lat), nts.tile250y(lat))
}


nts.bbox250 <- function(tile) {
  minlat <- tile[2] + 40
  wo <- nts.widthandoffset250(minlat)
  minlon <- -144 + wo[2] + (tile[1]*wo[1])
  maxlat <- minlat + 1
  maxlon <- minlon + wo[1]
  nts.makebbox(maxlat, maxlon, minlat, minlon)
}

nts.id250 <- function(tile250) {
  tileS <- nts.tileseriesfromtile250(tile250)
  seriesid <- nts.idseries(tileS)
  if(is.na(seriesid)) return(NA)
  
  seriesMinYTile <- tileS[2]*4
  yTileInSeries <- tile250[2]-seriesMinYTile
  
  mapsPerSeries <- nts.mapsperseries(tile250[2])
  seriesMinXTile <- tileS[1] * mapsPerSeries
  xTileInSeries <- tile250[1]-seriesMinXTile
  
  arealetter <- NULL
  if(tileS[2] >= 7) {
    arealetter <- nts.MAP_250K_N_OF_68[yTileInSeries+1, xTileInSeries+1]
  } else {
    arealetter <- nts.MAP_250K[yTileInSeries+1, xTileInSeries+1]
  }
  c(seriesid, arealetter)
}

nts.tile250byid <- function(ntsid) {
  if(length(ntsid) < 2) stop("Invalid NTS id in nts.tile250byid(): ", ntsid)
  seriesT <- nts.tileseriesbyid(ntsid[1])
  if(is.null(seriesT)) stop("Invalid NTS id in nts.tile250byid(): ", ntsid)
  if(seriesT[2]>=7) {
    map <- nts.MAP_250K_N_OF_68
  } else {
    map <- nts.MAP_250K
  }
  result <- nts.indexxy(ntsid[2], map)
  if(length(result)==0) stop("Invalid NTS id in nts.tile250byid() ", ntsid)
  tiley <- seriesT[2] * 4 + result[1] - 1
  mapsperseries <- nts.mapsperseries(tiley)
  tilex <- seriesT[1] * mapsperseries + result[2] - 1
  c(tilex, tiley)
}

nts.tile50byid <- function(ntsid) {
  tile250 <- nts.tile250byid(ntsid)
  if(length(ntsid) < 3) stop("Invalid NTS for 50k tile: ", ntsid)
  result <- nts.indexxy(ntsid[3], nts.MAP_50K)
  if(length(result)==0) stop("Invalid NTS: ", ntsid, " (", ntsid[3], ")")
  tiley <- tile250[2]*4+result[1]-1
  tilex <- tile250[1]*4+result[2]-1
  c(tilex, tiley)
}

nts.tile50y <- function(lat) {
  as.integer(floor((lat-40.0)/0.25))
}

nts.tile50x <- function(lon, lat) {
  tile250x <- nts.tile250x(lon, lat)
  wo <- nts.widthandoffset250(lat)
  londiff <- lon - (-144.0+wo[2]+(tile250x*wo[1]))
  plustilesx <- as.integer(floor(4.0*londiff/wo[1]))
  return(4*tile250x+plustilesx)
}

nts.tile50 <- function(lon, lat) {
  c(nts.tile50x(lon, lat), nts.tile50y(lat))
}

nts.bbox50 <- function(tile50) {
  minlat <- 40+tile50[2]*0.25
  wo <- nts.widthandoffset250(minlat)
  wd <- wo[1] / 4.0
  minlon <- -144.0 + wo[2] + tile50[1]*wd
  nts.makebbox(minlat+0.25, minlon+wd, minlat, minlon)
}


nts.id50 <- function(tile50) {
  tile250 <- as.integer(floor(tile50/4.0))
  id250 <- nts.id250(tile250)
  if(is.na(id250[1])) return(NA)
  
  mintiles <- tile250 * 4
  plustiles <- tile50 - mintiles
  sheet <- nts.MAP_50K[plustiles[2]+1, plustiles[1]+1]
  c(id250, sheet)
}


nts.bybboxgeneric <-function(bbox, tilefuncx, tilefuncy, idfunc) {
  minx <- max(bbox[1,1], -144.0)
  miny <- max(bbox[2,1], 40.0)
  maxx <- min(bbox[1,2], -48.0)
  maxy <- min(bbox[2,2], 88.0)
  if((maxx < minx) || (maxy < miny)) stop("Bounds provided may be outside the NTS grid")
  
  containsabove80 <- maxy > 80.0
  containsabove68 <- containsabove80 || (miny >= 68.0) || (maxy>68.0)
  containsbelow68 <- miny < 68.0
  containsbelow80 <- miny < 80.0
  
  idlist <- list()
  ld <- function(n=maxy, e=maxx, s=miny, w=minx) {
    mint <- c(tilefuncx(w, s), tilefuncy(s))
    maxt <- c(tilefuncx(e, n), tilefuncy(n))
    
    listind <- length(idlist)+1
    for(x in mint[1]:maxt[1]) {
      for(y in mint[2]:maxt[2]) {
        tileid <- idfunc(c(x,y))
        if(!is.na(tileid[1])) {
          idlist[[listind]] <<- tileid
        }
        listind <- listind+1
      }
    }
  }
  
  if(containsabove80) {
    if(containsbelow80) {
      ld(s=80)
    } else {
      ld()
      return(idlist)
    }
  }
  
  tempn <- maxy
  temps <- miny
  if(containsabove68) {
    if(containsabove80) {
      tempn <- 79.99
    }
    if(containsbelow68) {
      temps <- 68.0
    }
    ld(s=temps, n=tempn)
  }
  
  if(containsbelow68) {
    tempn <-maxy
    if(containsabove68) {
      tempn <- 67.99
    }
    ld(n=tempn)
  }
  idlist
}

# User friendly functions ----

#' A contstant denoting NTS Series scale (0)
#' @export
nts.SCALESERIES <- 0

#'  A constant denoting NTS Map Area (1:250k) scale (1)
#'  @export
nts.SCALE250K <- 1

#'  A constant denoting NTS Map Sheet (1:50k) scale (2)
#'  @export
nts.SCALE50K <- 2

#' Get NTS References by Bounding Box
#' 
#' Retreive a list of NTS references at a given scale by bounding box.
#' Bounding box is in the form returned by \code{sp::bbox()}. NTS References
#' all have a valid series component (e.g. "021"), but map area (e.g. "H")
#' and map sheet (e.g. "01") are not checked to make sure they exist.
#' 
#' @param bbox A bounding box of lat/lon values in the form returned by \code{sp::bbox()}.
#' @param atscale One of \code{nts.SCALESERIES}, \code{nts.SCALE250K}, or \code{nts.SCALE50K}
#' @return A \code{list} object containing zero or more NTS References.
#' @seealso \link{nts}
#' 
#' @export
nts.bybbox <- function(bbox, atscale) {
  if(atscale == nts.SCALESERIES) {
    tilexfunc <- nts.tileseriesx
    tileyfunc <- nts.tileseriesy
    idfunc <- nts.idseries
  } else if(atscale == nts.SCALE250K) {
    tilexfunc <- nts.tile250x
    tileyfunc <- nts.tile250y
    idfunc <- nts.id250
  } else if(atscale == nts.SCALE50K) {
    tilexfunc <- nts.tile50x
    tileyfunc <- nts.tile50y
    idfunc <- nts.id50
  } else {
    stop("Invalid value for atscale: ", atscale)
  }
  nts.bybboxgeneric(bbox, tilexfunc, tileyfunc, idfunc)
}

#' Get NTS Reference At A Location
#' 
#' Get NTS Reference(s) based on location at a given scale.
#' 
#' @param lat A scalar or vector of latitude values
#' @param lon A scalar or vector of longitude values of same length as \code{lat}
#' @param atscale One of \code{nts.SCALESERIES}, \code{nts.SCALE250K}, 
#' or \code{nts.SCALE50K}
#' @return A \code{list} of NTS References
#' 
#' @seealso \link{nts}
#' 
#' @export
nts.idat <- function(lat, lon, atscale) {
  if(length(lat) != length(lon)) stop("Length of lat and lon arguments must be the same in nts.idat")
  if(atscale == nts.SCALESERIES) {
    tilefunc <- nts.tileseries
    idfunc <- nts.idseries
  } else if(atscale == nts.SCALE250K) {
    tilefunc <- nts.tile250
    idfunc <- nts.id250
  } else if(atscale == nts.SCALE50K) {
    tilefunc <- nts.tile50
    idfunc <- nts.id50
  } else {
    stop("Invalid value for atscale: ", atscale)
  }
  out <- list()
  for(i in 1:length(lat)) {
    out[[i]] <- idfunc(tilefunc(lon[i],lat[i]))
  }

  out
}

#' Get Bounding Box of an NTS Reference
#' 
#' Calculates the bounding box in latitude/longitude described by
#' a particular NTS Reference.
#' 
#' @param ntsid One or more NTS References as generated by \code{nts()}
#' @return A bbox like that returned by \code{sp::bbox()}, or a list of
#' such objects.
#' @examples nts.bbox(nts('21h'))
#' 
#' @export
nts.bbox <- function(ntsid) {
  if(class(ntsid)=="list") {
    out <- list()
    for(i in 1:length(ntsid)) {
      out[[i]] <- nts.bbox(ntsid[[i]])
    }
    out
  } else {
    if(length(ntsid)>=3) {
      fun <- nts.bbox50
      tilef <- nts.tile50byid
    } else if(length(ntsid)==2) {
      fun <- nts.bbox250
      tilef <- nts.tile250byid
    } else if(length(ntsid)==1) {
      fun <- nts.bboxseries
      tilef <- nts.tileseriesbyid
    } else {
      stop("Invalid NTS: ", ntsid)
    }
    tile <- tilef(ntsid)
    fun(tile)
  }
}
 
