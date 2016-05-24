# Copyright (c) 2003-4 by Barry Rowlingson and Roger Bivand

#.valid.CRSobj <- function(object) {
#	if (exists("is.R") && is.function(is.R) && is.R()) {
#		projargNA <- is.na(object@projargs)
#	} else {
#		projargNA <- is.na(as.numeric(object@projargs))
#	}
#	if (!projargNA) {
#		res <- .Call("checkCRSArgs", object@projargs, 
#				PACKAGE="rgdal")
#	} else res <- list(TRUE, as.character(NA))
#	if (!res[[1]]) {
#	    	return(res[[2]])
#	} else {
#		return(res[[1]])
#	}
#
#}
#
#setValidity("CRS", .valid.CRSobj)


"CRSargs" <- function(object) {
	if (!is(object, "CRS")) stop("not a CRS object")

	if (!is.na(object@projargs)) {
		res <- (checkCRSArgs(object@projargs)[[2]])
                res <- paste(unique(unlist(strsplit(res, " "))), 
			collapse=" ")
                return(sub("^\\s+", "", res))
	} else return(as.character(NA))
}

checkCRSArgs <- function(uprojargs) {
# RSB 2015-05-21
# fix for omission of proj_defs.dat in PROJ.4 4.9.1
  if (!get("has_proj_def.dat", envir=.RGDAL_CACHE)) {
      message("NOTE: rgdal::checkCRSArgs: no proj_defs.dat in PROJ.4 shared files")
      uprojargs <- proj_def_bug_fix(uprojargs)
  }
  res <- .Call("checkCRSArgs", uprojargs, PACKAGE="rgdal")
  res[[2]] <- sub("^\\s+", "", res[[2]])
# fix for pj_get_def() +no_uoff/+no_off bug
  no_uoff <- length(grep("+no_uoff", uprojargs, fixed=TRUE) > 0)
  no_off <- length(grep("+no_off", uprojargs, fixed=TRUE) > 0)
  if (no_uoff) {
      if( length(grep("+no_uoff", res[[2]], fixed=TRUE)) == 0) 
          res[[2]] <- sub("+no_defs", "+no_uoff +no_defs", res[[2]], fixed=TRUE)
  }
  if (no_off) {
      if (length(grep("+no_off", res[[2]], fixed=TRUE)) == 0)
          res[[2]] <- sub("+no_defs", "+no_off +no_defs", res[[2]], fixed=TRUE)
  }
  res
}

proj_def_bug_fix <- function(uprojargs) {
    if (length(grep("no_defs", uprojargs)) == 0L && 
# corrected 20150904
        length(grep("init", uprojargs)) == 0L) {
        if (length(grep("ellps", uprojargs)) == 0L && 
# corrected 20150905
            length(grep("datum", uprojargs)) == 0L) {
            tags <- sapply(strsplit(strsplit("+proj=longlat +no_defs",
                "\\+")[[1]], "="), "[", 1)
# based on proj/src/pj_init.c lines 191-197
            if (!any(c("datum", "ellps", "a", "b", "rf", "f") %in% tags)) {
                uprojargs <- paste(uprojargs, "+ellps=WGS84", sep=" ")
           }
         }
    }
    uprojargs
}
