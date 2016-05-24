"timeSeries" <- 
function(data, positions., units., from = timeCalendar(d = 1, m = 1, y = 1960),
	by = "days", k.by = 1, align.by = FALSE, week.align = NULL)
{
  ## function to create a timeSeries
  if(missing(positions.) && missing(data) && missing(units.) && missing(from)
     && missing(by)) return(new("timeSeries"))
  if(!missing(positions.) && (length(positions.) != nrow(as.data.frame(data))))
    stop("Positions and data lengths do not agree")
  ret <- new("timeSeries")
  ret@data <- asSeriesData(data)
  if(missing(positions.)) {
    len <- numRows(ret@data)
    from <- as(from, "timeDate")
    ret@positions <- timeSequence(from = from, length.out = len, by = by,
                                  k.by = k.by, align.by = align.by,
                                  week.align = week.align,
                                  format = from@format, zone = from@time.zone
                                  )
  } else {
    if(!is(positions., "positionsCalendar"))
      positions. <- as(positions., "timeDate")
    ret@positions <- positions.
  }
  if(!missing(units.))
    ret@units <- as(units., "character")
  ret
}

"positions" <- 
function(object)
{
	# return the positions of an ordered data object
	object@positions
}

"positions<-" <- 
function(object, value)
{
	# replace the positions of an ordered data object
	if(length(value) != numRows(object@data)) stop(
			"Positions and data lengths do not agree")
	if(!is(value, "positions") && !is(value, class(object@positions)))
		stop("object is not a valid positions object")
	object@positions <- value
	object
}

setMethod("seriesLag", signature(X="series"), 
          function(X, k = 1, trim = FALSE, pad = NA)
          {
            x <- seriesData(X)
            pos = positions(X)
            if(is.null(dim(x))) {
              isVec = TRUE
            }
            else {
              isVec = FALSE
            }
            x = as.matrix(x)
            tmp = dim(x)
            n = tmp[1]
            p = tmp[2]
            lk = length(k)
            if(lk == 0)
              stop("The optional argument k must have positive length.")
            tmp = matrix(pad, n, p * lk)
            storage.mode(tmp) = "double"
            ans = .C("series_lag",
              as.double(x),
              as.integer(n),
              as.integer(p),
              as.integer(k),
              as.integer(lk),
              mat = tmp,
              NAOK = TRUE)$mat
            dnames = colIds(x)
            if(length(dnames) == 0) {
              if(p > 1)
                dnames = paste("V", 1:p, sep = "")
              else dnames = ""
            }
            k.lab = k
            if(p == 1 && dnames == "") {
              tmp = rep("lag", lk)
              if(any(k < 0)) {
                tmp[k < 0] = "lead"
                k.lab[k < 0] =  - k[k < 0]
              }
              dnames = paste(tmp, k.lab, sep = "")
            }
            else {
              tmp = rep("lag", lk)
              if(any(k < 0)) {
                tmp[k < 0] = "lead"
                k.lab[k < 0] =  - k[k < 0]
              }
              dnames = outer(dnames, paste(tmp, k.lab, sep = ""), FUN = paste,
                sep = ".")
            }
            dimnames(ans) = list(NULL, as.vector(dnames))
            if(trim) {
              if(any(k < 0)) {
                kn = k[k < 0]
                tmp = (n + min(kn) + 1):n
                kp = k[k >= 0]
                if(length(kp) > 0) {
                  tmp = c(tmp, 1:max(kp))
                }
              }
              else {
                tmp = 1:max(k)
              }
              ans = ans[ - tmp,  , drop = FALSE]
              if(isVec && lk == 1) {
                ans = as.vector(ans)
              }
              if(inherits(X, "timeSeries")) {
                ans = timeSeries(ans, positions. = pos[ - tmp])
              } else if (inherits(X, "signalSeries")) {
                ans = signalSeries(ans, positions. = pos[ - tmp])
              }
            }
            else {
              if(isVec && lk == 1)
                ans = as.vector(ans)
              if(inherits(X, "timeSeries")) {
                ans = timeSeries(ans, positions. = pos)
              } else if (inherits(X, "signalSeries")) {
                ans = signalSeries(ans, positions. = pos)
              } 
            }
            ans
          })
