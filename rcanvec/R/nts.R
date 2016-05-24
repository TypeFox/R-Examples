#NTS tools (translated from Java)

#nts format is a vector of character strings, e.g. c("021", "h", "01")
#requires source of ntstiles.R

.isanum <- function(x) {
  suppressWarnings(!is.na(as.numeric(x)))
}

.validseries <- function(x) {
  .isanum(x)
}

.validarea <- function(x) {
  !.isanum(x) && nchar(x) == 1
}

.validsheet <- function(x) {
  if(.isanum(x)) {
    if(1 <= as.numeric(x) && as.numeric(x) <= 16 ) {
      return(TRUE)
    }
  }
  FALSE
}

.addzeroes <- function(x, len) {
  if(nchar(x) >= len) {
    return(x)
  }
  out <- x
  for(i in 1:(len-nchar(x))) {
    out <- paste0("0", out)
  }
  out
}

.parsentsstring <- function(st) {
  st <- as.character(st)
  series <- NULL
  area <- NULL
  sheet <- NULL
  
  for(i in 1:nchar(st)) {
    subs <- substr(st, 1, i)
    char <- substr(st, i, i)
    if(!.isanum(subs) && i==1) {
      stop("Invalid NTS: ", st)
    } else if(i>1 && !.isanum(char)) {
      #substring is no longer a number, so area letter must be current char
      series <- substr(st, 1, i-1)
      if(!.validseries(series)) {
        stop("Invalid NTS: ", st, " (", series, ")")
      }
      break
    } 
  }
  
  if(is.null(series)) {
    #whole argument may be a series
    if(.validseries(st)) {
      series <- st
    } else {
      stop("Invalid NTS: ", st, " (", st, ")")
    }
  }
  
  #if series was not the whole argument, try to find the area letter
  if(nchar(series) < nchar(st)) {
    area <- substr(st, nchar(series)+1, nchar(series)+1)
    if(!.validarea(area)) {
      stop("Invalid NTS Sheet: ", st, " (", area, ")")
    }
    
    if((nchar(series) + nchar(area)) < nchar(st)) {
      sheet <- substr(st, nchar(series)+nchar(area)+1, nchar(st))
      if(!.validsheet(sheet)) {
        stop("Invalid NTS Sheet: ", st, " (", sheet, ")")
      }
    }
  }
  
  
  #if series and area is not the whole argument, sheet is the remainder
  
  
  series <- .addzeroes(series, 3)
  if(is.null(area)) {
    return(series)
  } else {
    area <- toupper(area)
    if(is.null(sheet)) {
      return(c(series, area))
    } else {
      sheet <- .addzeroes(sheet, 2)
      return(c(series, area, sheet))
    }
  }
  
  stop("Error in .ntsstring(), reached end of function by error")
}

#' Generate NTS References
#' 
#' Generate one or more NTS references based on arguments provided.
#' 
#' @param ... An arbitrary number of strings in the form 21h1 to be parsed
#' @param lat A vector of latitude values
#' @param lon A vector of longitude values
#' @param bbox A bounding box in the form returned by sp::bbox()
#' @param atscale An integer value describing scale, either 1 (250k) or 2 (50k)
#' @return one or more NTS references in the form c("021", "H", "01")
#' @examples
#' nts('21h')
#' nts('21h1')
#' nts('21h1', '21a16', '021A15')
#' nts(lat=45.2, lon=-64.32)
#' nts(lat=c(45.2, 46.2), lon=c(-64.32, -64.81))
#' 
#' library(prettymapr)
#' nts(bbox=makebbox(45.125, -64.25, 44.875, -64.75))
#' 
#' @export
#' 
nts <- function(..., lat=NULL, lon=NULL, bbox=NULL, atscale=nts.SCALE50K) {
  ntsstring <- list(...)
  
  if(!is.null(bbox)) {
    if(length(ntsstring)>0) stop("bbox and nts strings both specified to nts()")
    if(length(lat)>0 || length(lon)>0) stop("bbox and lat/lon pairs both specified to nts()")
    out <- nts.bybbox(bbox, atscale)
  } else if(length(lat)>0 || length(lon)>0) {
    if(length(ntsstring)>0) stop("lat/lon lists and nts strings both specified to nts()")
    if(length(lat) != length(lon)) stop("length of lat and lon lists must be the same")
    out <- nts.idat(lat, lon, atscale)
  } else if(length(ntsstring)>0) {
    out <- list()
    for(i in 1:length(ntsstring)) {
      #in case an argument is a characters vector like that put out by ntsstring
      arg <- ntsstring[[i]]
      for(j in 1:length(arg)) {
        out[[length(out)+1]] <- .parsentsstring(arg[j])
      }
    }
    
  } else {
    stop("No arguments passed to nts()")
  }
  
  if(length(out)==1) {
    out[[1]]
  } else {
    out
  }
}

#' Generate NTS Reference Strings
#' 
#' Generate NTS strings from NTS references or other arguments to generate
#' such a list from nts() based on arguments provided.
#' 
#' @param ntsid An arbitrary number of NTS references in the form c("021", "H", "01") 
#' or \code{list(c("021", "H", "01"), c("021", "A", "16"))} to be converted
#' @param lat A vector of latitude values
#' @param lon A vector of longitude values
#' @param bbox A bounding box in the form returned by sp::bbox()
#' @param atscale An integer value describing scale, either 1 (250k) or 2 (50k)
#' @return a character vector of NTS reference strings in the form "021H01".
#' @examples
#' ntsstring(c("021", "H", "01"))
#' 
#' library(prettymapr)
#' ntsstring(bbox=makebbox(45.125, -64.25, 44.875, -64.75))
#' 
#' @export
#' 
ntsstring <- function(ntsid=NULL, lat=NULL, lon=NULL, bbox=NULL, atscale=nts.SCALE50K) {
  if(!is.null(lat) || !is.null(lon) || !is.null(bbox)) {
    if(!is.null(ntsid)) stop("NTS ids specified alongside arguments to be passed to nts()")
    return(ntsstring(nts(lat=lat, lon=lon, bbox=bbox, atscale=atscale)))
  } else if(!is.null(ntsid)) {
    if(class(ntsid)=="list") {
      out <- rep("", length(ntsid))
      for(i in 1:length(ntsid)) {
        out[i] <- paste0(ntsid[[i]], collapse="")
      }
      out
    } else {
      paste0(ntsid, collapse="")
    }
  } else {
    stop("No arguments passed to ntsstring()")
  }

}

