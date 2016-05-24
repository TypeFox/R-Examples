
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  




earth_coordinate_names<- function(names) {
  ## Earth coordinates + possibly radius
  n <- substr(tolower(names), 1, 6)
  nc <- nchar(n)
  lon <- lat <- logical(length(n))
  for (i in 1:length(n)) {
    lon[i] <- substr("longit", 1, nc[i]) == n[i]
    lat[i] <- substr("latitu", 1, nc[i]) == n[i]
  }
  lonORlat <- lon | lat  
  earth <- all(nc[lonORlat] >= 2) && sum(lon==1) && sum(lat == 1)
  
  return(if (length(names)==2 | !earth) earth else         
         if ( (lo <- which(lon)) < (la <- which(lat)))
         which(lonORlat[]))
}

cartesian_coordinate_names <- function(names) {
  n <- substr(tolower(names), 1, 1)
  coords <- c("T", "x", "y", "z")
  Txyz <- outer(n, coords, "==")
  cs <- colSums(Txyz)
  if (any(cs > 1) || sum(cs[1:2]) == 0 || any(diff(cs[-1]) > 0))
    return (integer(0))
  Txyz <- Txyz[, c(2:4, 1), drop=FALSE]
  ord <- apply(Txyz, 2, function(x) which(x > 0))
  ord <- order(unlist(ord))
  rs <- which(rowSums(Txyz) > 0)
  return(rs[ord])
}



general_coordinate_names <- function(names) {
  n <- substr(tolower(names), 1, 5)
  return(which(n == "coord"))
}



data.columns <- function(data, xdim=0, force=FALSE, halt=TRUE) {
  #  Print("data.col", data)
   if (length(xdim) == 0) xdim <- 0
  if (xdim>0 && xdim >= ncol(data)) stop("not enough columns in 'data'.")
  RFopt <- RFoptions()
  info <- RFoptions()$coords
  cn <- colnames(data)

  if (all(is.na(info$varnames))) {
    if (all(is.na(info$coordnames))) {
      if (is.null(cn)) {
        if (force) return(list(data=(xdim+1):ncol(data), x=1:xdim))
        if (halt)
          stop('colnames of data argument must contain "data" or "variable"')
        else return(NULL);
      }
      is.data <- (tolower(substr(cn, 1, 4)) == "data" |
                  tolower(substr(cn, 1, 4)) == "value" |
                  tolower(substr(cn, 1, 8)) == "variable")
      if (!any(is.data)) {
        if (force) return(list(data=(xdim+1):ncol(data), x=1:xdim))
        if (halt) stop('no colname starts with "data" or "variable"')
        else return(NULL);
      }
      is.data <- which(is.data)
      if (is.data[1] > 1) is.data <- is.data[1] : ncol(data)# coord am Anfang
      ##     dann wird Rest alles als data angenommen, egal welcher Name
    }
  } else {
    if (is.numeric(info$varnames)) {
      is.data <- rep(info$varnames, length.out=2)
      if (is.na(is.data[1])) is.data[1] <- 1
      if (is.na(is.data[2])) is.data[2] <- ncol(data)
      is.data <- is.data[1] : is.data[2]
    } else {
      if (RFopt$general$vdim_close_together)
        stop("'vdim_close_together' must be FALSE")
      l <- list()
      vdim <- length(info$varnames)
      for (v in 1:vdim)
        l[[v]] <-
          substring(cn, 1, nchar(info$varnames[v])) == info$varnames[v]
      repet <- sapply(l, sum)
      if (repet[1] == 0) stop("data names could not be detected") 
      if (any(repet != repet[1]))
        stop("detected repetitions are not equal for all components")
      m <- matrix(unlist(l), ncol=vdim)
      if (any(rowSums(m) > 1))
        stop("names of multivariate components not unique")
      is.data <- as.vector(t(apply(m, 2, which)))
    }
  }

  if (all(is.na(info$coordnames))) {
    is.x <- (1:ncol(data))[-is.data]
    if (xdim > 0) {
      if (length(is.x) < xdim)
        stop("not enough columns for coordinates found ")
      if (xdim < length(is.x) &&
          RFopt$general$printLevel >= PL_SUBIMPORTANT)
        message("column(s) '", paste(is.x[-1:-xdim], collapse=","),
                "' not used.\n")
      is.x <- is.x[1:xdim]
    }
  } else {
    if (is.numeric(info$coordnames)) {
      is.x <- rep(info$coordnames, length.out=2)
      if (is.na(is.x[1])) is.x[1] <- 1
      if (is.na(is.x[2])) is.x[2] <- ncol(data)
      is.x <- is.x[1] : is.x[2]
    } else {
      l <- list()
      len <- length(info$coordnames)
      for (i in 1:len)
        l[[v]] <-
          substring(cn, 1, nchar(info$coordnames[v])) == info$coordnames[v]
      is.x <- unlist(l)
      if (xdim > 0 && xdim != length(l))
        stop("expected dimension of coordinates does not match the found coordinates")
    }
     
    if (all(is.na(info$varnames))) {
      is.data <-  (1:ncol(data))[-is.x]
      if (length(is.data) == 0) stop("no columns for data found")
    } else {
     if (any(is.x %in% is.data))
       stop("column names and data names overlap.")
     if (length(is.x) + length(is.data) < ncol(data) &&
         RFopt$general$printLevel >= PL_SUBIMPORTANT)
       message("column(s) '",
               paste(1:ncol[c(-is.x, -is.data)], collapse=","),
               "' not used.\n")
    }
  }
  return(list(data=is.data, x=is.x) )
}

    
GetDataNames <- function(model, coords=NULL, locinfo) {#, data=NULL) {
#  if (length(data) > 0) {
#    cd <- try(Check Data(model=model, given=coords, data=data))      
#    if (class(cd) != "try-error")
#      return(list(coordnames=cd$coordnames, varnames=cd$varnames))
#  }

 #  Print(missing(model), missing(coords), missing(locinfo))
 # Print(model);  str(coords);  Print(coords);  Print(locinfo)

  varnames <- extractVarNames(model)
  coordnames <-
    if (!is.null(coords) && is.matrix(coords)) colnames(coords) else NULL
                       # gets response part of model, if model is a formula

  Zeit <- locinfo$Zeit

#  Print(locinfo)
  
  tsdim <- locinfo$spatialdim + locinfo$Zeit

#  Print(tsdim)
  
  if (is.null(coordnames)) {
    system <- RFoptions()$coords$coord_system
    if (system == "earth") {
      coordnames <- if (tsdim == 4) ZF_EARTHCOORD_NAMES
      else if (tsdim == 2) ZF_EARTHCOORD_NAMES[1:2]
      else if (tsdim == 3) c(ZF_EARTHCOORD_NAMES[1:2], "HeightOrTime")
    } else if (system == "cartesian" && tsdim <= 4) {
      coords <- ZF_CARTCOORD_NAMES[1:tsdim]
      if (Zeit) coordnames[tsdim] <- ZF_CARTCOORD_NAMES[3 + 1]
    } else {
      coordnames <- paste(ZF_GENERAL_COORD_NAME[1], 1:tsdim, sep="")
      if (Zeit) coordnames[tsdim] <- ZF_GENERAL_COORD_NAME[2]
    }
  }

  return(list(coordnames=coordnames, varnames=varnames))
}



search.model.name <- function(cov, name, level) {
  if (length(name) == 0 || length(cov) ==0) return(cov);
  if (!is.na(pmatch(name[1], cov))) return(search.model.name(cov, name[-1], 1))

  for (i in 1:length(cov$submodels)) {
    found <- search.model.name(cov$submodels[[i]], name, 1)
    if (!is.null(found)) return(found)      
  }
  found <- search.model.name(cov$internal, name, 1)
  if (!is.null(found)) return(found)
  if (level == 0) stop("model name not found")
  return(NULL)
}



