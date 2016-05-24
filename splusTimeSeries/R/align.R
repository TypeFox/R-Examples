"align" <- 
function(x, pos, how = "NA", error.how = "NA", localzone = FALSE, matchtol = 0,
	by, k.by = 1, week.align = NULL, holidays = timeDate())
{
  ## Align series object x to new positions
  if((how != "NA") && (how != "drop") && (how != "nearest") &&
     (how != "before") && (how != "after") && (how != "interp"))
    stop(paste("unknown value \"", how, "\"for how argument", sep = ""))
  if((error.how != "NA") && (error.how != "drop") && (error.how != "nearest"))
    stop(paste("unknown value \"", error.how, 
               "\"for error.how argument", sep = ""))
  if(matchtol < 0)
    stop("matchtol must be >= 0")
  ## construct positions if missing
  origpos <- positions(x)
  if(missing(pos)) {
    rng <- range(origpos)
    if(is(by, "character"))
      pos <- timeSeq(from = rng[1], to = rng[2], by = by,
                     k.by = k.by, align.by = TRUE, extend = TRUE, 
                     week.align = week.align, holidays = holidays)
    else pos <- timeSeq(from = rng[1], to = rng[2], by = by)
  }
  else {
    if(length(pos) > 1 && any(diff(pos) < 0)) {
      warning("pos is not sorted, results will be sorted")
      pos <- sort(pos)
    }
  }
  if(identical(origpos, pos))
    return(x)
  origlen <- length(origpos)
  newlen <- length(pos)
  newpos <- pos
  if(!newlen || ( !origlen && error.how == "drop" )) {
    if( !origlen ) return( x )
    else return(x[0])
  }
  if( !origlen ) {
    ## original sequence was zero length, but new positions aren't,
    ## and the person didn't want to drop missing rows
    numRows( x@data ) <- newlen
    x@positions <- pos
    return( x )
  }
  if(localzone) {
    ## efficient way to convert so comparison is done in local zone is
    ## to go to character (w/out time zone) and read back as default zone
    ## (timeZoneConvert goes the other direction)
    origpos@format <- "%m %d %Y %H %M %S %N"
    newpos@format <- "%m %d %Y %H %M %S %N"
    origpos <- timeDate(as(origpos, "character"), in.format = 
			"%m %d %Y %H %M %S %N")
    newpos <- timeDate(as(newpos, "character"), in.format = 
                       "%m %d %Y %H %M %S %N")
  }
  ## call appropriate align function from splusTimeDate to find out alignment
  ## returns a list of (indexes for NA, indexes to drop, subscript or
  ## interpolation weights)
  if(is(pos, "positionsCalendar")) {
    how.align <- splusTimeDate:::timealign(origpos, newpos, how, error.how,
        matchtol)
  }
  else {
    ##pos = as(pos, "numericSequence")
    how.align <- splusTimeDate:::numalign(origpos, newpos, how, error.how,
        matchtol)
  }
  keeprows <- !how.align[[2]]
  narows <- how.align[[1]] & keeprows
  needrows <- keeprows & !narows
  newdata <- x@data
  olddata <- x@data
  ## get the right length on data before subscripting
  numRows(newdata) <- newlen
  if(how == "interp") {
    ## be careful -- S doesn't seem to allow direct replacement
    ## of integers with numerics here!
    newd <- how.align[[3]][needrows] *
      subscript2d(olddata, how.align[[4]][needrows],  ) 
    ## don't add in places where weight is zero, to avoid NA issues
    otherwt <- how.align[[5]][needrows]
    otherdat <- subscript2d( olddata, how.align[[6]][needrows], )
    whereadd <- otherwt != 0
    if( any( whereadd )) {
      subscript2d(newd, whereadd, ) <- subscript2d( newd, whereadd, ) +
        otherwt[ whereadd ] * subscript2d( otherdat, whereadd, )
    }
    if((class(newd) == "numeric") && (class(newdata) == "integer")
       )
      newdata <- as(newdata, "numeric")
    subscript2d(newdata, needrows,  ) <- newd
  }
  else subscript2d(newdata, needrows,  ) <-
    subscript2d(olddata, how.align[[3]][needrows],
			)
  ## 
  if(any(narows)) subscript2d(newdata, narows,  ) <- rep(NA, sum(narows))
  if(!all(keeprows)) {
    x@data <- asSeriesData(subscript2d(newdata, keeprows,  ))
    x@positions <- pos[keeprows]
  }
  else {
    x@data <- asSeriesData(newdata)
    x@positions <- pos
  }
  x
}

"seriesMerge" <- 
function(x1, x2, ..., pos = positions(x1), how, error.how, localzone = FALSE, 
	matchtol = 0, suffixes)
{
  ## pos can be "union" or a positions object.
  ## if it is "union", then how and error.how should not be set to "drop", 
  ## and default to "NA".  Otherwise, they default to "drop" and the
  ## default args make intersection positions of all the inputs.
  otherarglist <- list(x2, ...)
  ndots <- length(otherarglist) - 1
  allarglist <- list(x1, x2, ...)
  ## figure out if we want union of all positions or passed-in positions
  posunion <- identical(pos, "union")
  if(posunion) {
    if(missing(error.how))
      error.how <- "NA"
    if(missing(how))
      how <- "NA"
    if(error.how == "drop") {
      warning(paste(
                    "error.how=\"drop\" does not make sense with",
                    "positions=\"union\", setting to NA"))
      error.how <- "NA"
    }
    if(how == "drop") {
      warning(paste("how=\"drop\" does not make sense with",
                    "positions=\"union\", setting to NA"))
      how <- "NA"
    }
  }
  else {
    if(missing(error.how))
      error.how <- "drop"
    if(missing(how))
      how <- "drop"
  }
  ## figure out ncol for each arg, will need below
  ncollist <- sapply(allarglist, "numCols")
  ## get suffixes default, will need below
  if(missing(suffixes)) suffixes <- paste(".", 1:(ndots + 2), sep = "")
  if((length(suffixes) != (ndots + 2)) || any(duplicated(suffixes)))
    stop(paste("Must be", ndots + 2, " distinct suffixes"))
  ## define a function to merge two args; names may be bogus on return
  merge2 <- function(x, y, how, error.how, localzone, matchtol)
    {
      y <- align(y, positions(x), how, error.how, localzone, matchtol
                 )
      x <- align(x, positions(y), how, error.how, localzone, matchtol
                 )
      x@data <- asSeriesData(cbind(x@data, y@data))
      x@units <- c(x@units, y@units)
      x
    }
  ## merge the positions, if we're supposed to
  if(posunion) {
    posl <- lapply(allarglist, "positions")
    posl[["localzone"]] <- localzone
    posl[["matchtol"]] <- matchtol
    pos <- do.call("unionPositions", posl)
  }
  ## now align the first arg to the calculated or input positions
  x1 <- align(x1, pos, how, error.how, localzone, matchtol)
  ids <- colIds(x1)
  if(length(ids) == 0) {
    ids <- seq(len=numCols(x1))
  }
  nms <- ids
  ## and then merge in everything else
  for(arg in otherarglist) {
    x1 <- merge2(x1, arg, how, error.how, localzone, matchtol)
    ids <- colIds(arg)
    if(length(ids) == 0) {
      ids <- seq(len=numCols(arg))
    }
    nms <- c(nms, ids)
  }
  ## fix the column names
  dups <- unique(nms[duplicated(nms)])
  if(length(dups)) {
    ## make a vector of which arg each name would have come from
    ncolvec <- unlist(lapply(1:length(ncollist),
                             function(i, ncollist)
                             rep(i, ncollist[[i]]), ncollist))
    for(d in dups)
      nms[nms == d] <- paste(d, suffixes[ncolvec[nms == d]],
            sep = "")
  }
  ## may have ended up with null matrix here
  if(is.null(x1@data)) rowIds(x1@data) <- character(0)
  if(length(nms) && ncol( x1@data ))
    colIds(x1@data) <- nms
  x1
}

