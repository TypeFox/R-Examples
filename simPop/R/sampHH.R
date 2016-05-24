#' Sample households from given microdata.
#' 
#' The function samples households from microdata containing personal and
#' household information.
#' 
#' households are drawn from the data and new ID's are generated for the new
#' households.
#' 
#' @name sampHH
#' @param pop data frame containing households and persons
#' @param sizefactor factor of how many times the initial population should be
#' resampled
#' @param hid string specifying the name of the household-id variable in the
#' data.
#' @param strata can be used to sample within strata.
#' @param hsize string specifying the name of the household size variable in
#' the data.
#' @return the data frame of new households.
#' @export
#' @author Bernhard Meindl, Matthias Templ and Johannes Gussenbauer
#' @keywords manip
#' @examples
#' data(eusilcP)
#' pop <- eusilcP
#' colnames(pop)[3] <- "hhsize"
#' 
#' system.time(x1 <- sampHH(pop, strata="region", hsize="hhsize"))
#' dim(x1)
#' \dontrun{
#' ## approx. 10 second computation time ...
#' system.time(x1 <- sampHH(pop, sizefactor=4, strata="region", hsize="hhsize"))
#' dim(x1)
#' system.time(x2 <- sampHH(pop, strata=NULL, hsize="hhsize"))
#' 
#' pop <- pop[,-which(colnames(pop)=="hhsize")]
#' system.time(y1 <- sampHH(pop, strata="region", hsize=NULL))
#' system.time(y2 <- sampHH(pop, strata=NULL, hsize=NULL))
#' }
sampHH <- function(pop, sizefactor=1, hid="hid", strata="region", hsize=NULL) {
  x <- nr <- id <- NULL
  pop <- as.data.table(pop)
  setkeyv(pop, hid)

  if ( is.null(hsize) ) {
    pop$x <- 1
    pop[,hsize:=list(hsize=sum(x)), by=hid][,x:=NULL]
    tmphsize <- "hsize"
  } else {
    if ( is.na(match(hsize, colnames(pop))) ) {
      stop("provide a valid variable containing household sizes!\n")
    }
    tmphsize <- hsize
  }
  rm(hsize)
  hsize.ind <- match(tmphsize, colnames(pop))

  if ( !is.null(strata) && is.na(match(strata, colnames(pop))) ) {
    stop("provide a valid variable containing stratas!\n")
  }

  if ( !is.null(strata) ) {
    frame <- pop[,c(hid,strata), with=FALSE]
    setkeyv(frame, hid)
    frame <- unique(frame)
    s <- split(frame, frame[,match(strata, colnames(frame)),with=FALSE])
    ids <- lapply(s, function(x) {
      res <- data.table(hid=sample(x[[hid]], round(nrow(x)*sizefactor), replace=TRUE))
      setnames(res, hid)
      res
    })

    xx <- lapply(ids, function(x) {
      merge(x, pop, allow.cartesian=TRUE, by=hid, suffixes = c("_x", "_y"))
    })
    out <- rbindlist(xx)
    ids <- rbindlist(ids)
  } else {
    ids <- data.table(hid=sample(pop[,hid], round(nrow(pop)*sizefactor), replace=TRUE))
    out <- merge(ids, pop, allow.cartesian=TRUE, by=hid, suffixes = c("_x", "_y"))
  }

  # recode hsize-variable to numeric
  out[[hsize.ind]] <- as.numeric(out[[tmphsize]])
  setkeyv(out, hid)

  # nicht optimal, geht sicher einfacher
  setkeyv(ids, hid)
  ids$x <- 1
  ids <- ids[,.N,by=hid]

  ids$hsize <- out[,list(x=max(get(tmphsize))), by=hid]$x
  ids <- data.table(id=unlist(sapply(1:nrow(ids), function(x) {
    rep(paste(ids[[hid]][x], 1:ids[["N"]][x], sep="."), each=ids[["hsize"]][x])
  })))
  ids$x <- 1
  ids <- ids[,list(sum(x)),by=id]
  out[[hid]] <- rep(1:length(unique(ids$id)), times=ids$V1)
  out$hid <- out[[hid]] + max(pop[[hid]])
  invisible(out)
}

