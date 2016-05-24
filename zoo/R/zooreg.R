zooreg <- function(data, start = 1, end = numeric(), frequency = 1, 
  deltat = 1, ts.eps = getOption("ts.eps"), order.by = NULL)
{
    ## choose frequency/deltat
    if (missing(frequency)) frequency <- 1/deltat
    	else if(missing(deltat)) deltat <- 1/frequency
    if (frequency > 1 && abs(frequency - round(frequency)) < ts.eps)
    	frequency <- round(frequency)

    ## check data and choose default
    if (missing(data) || is.null(data)) data <- NA
    if(!(is.vector(data) || is.factor(data) || is.atomic(data) || is.matrix(data) || is.data.frame(data)))
      stop(paste(dQuote("data"), ": attempt to define invalid zoo object"))
    if(is.matrix(data) || is.data.frame(data)) data <- as.matrix(data)

    ## if no index (i.e., order.by) is specified: behave as ts()
    ## else: behave as zoo()
    if (is.null(order.by)) {
	if(!any(c(is.vector(data), is.factor(data), is.atomic(data), is.matrix(data), is.data.frame(data))))
  	    stop(paste(dQuote("data"), ": attempt to define invalid zoo object"))
	ndata <- NROW(data)        

        ## convenience function
        numORint <- function(x) identical(class(x), "numeric") | identical(class(x), "integer")

        ## choose start/end
        if (length(start) > 1) start <- start[1] + (start[2] - 1)/frequency
        if (length(end) > 1) end <- end[1] + (end[2] - 1)/frequency
        if (missing(end)) {
	  ostart <- start
	  oend <- NULL
	  start <- as.numeric(start)	  
	  end <- start + (ndata - 1)/frequency
	} else if(missing(start)) {
	  ostart <- NULL
	  oend <- end
	  end <- as.numeric(end)
	  start <- end - (ndata - 1)/frequency
	} else{
	  ostart <- start
	  oend <- NULL
	  start <- as.numeric(start)
	  end <- as.numeric(end)
	}
        if (start > end) stop("start cannot be after end")

        ## check whether lengths of data and index match
	# wrong lengths sometimes: order.by <- seq(start, end, by = deltat)
	order.by <- start + seq(0, length.out = ndata) * deltat
	if(identical(all.equal(start*frequency, round(start*frequency)), TRUE)) {
	  order.by <- floor(frequency*order.by + .0001)/frequency
        }
	
	## support also non-numeric indexes
	if(!is.null(ostart) && !numORint(ostart))
	  order.by <- ostart + (order.by - start)
	if(!is.null(oend) && !numORint(oend))
	  order.by <- oend + (order.by - end)
	
	nobs <- length(order.by)
        ## nobs <- floor((end - start) * frequency + 1.01)
        if (nobs != ndata) {
	  if(is.vector(data)) data <- rep(data, length.out = nobs)
	  else if(is.factor(data)) data <- factor(rep(as.character(data), length.out = nobs), labels = levels(data))
	  else if(is.matrix(data) || is.data.frame(data)) data <- data[rep(1:ndata, length.out = nobs), , drop = FALSE]
        }
 
        attr(data, "oclass") <- attr(data, "class")
        attr(data, "index") <- order.by
        attr(data, "frequency") <- frequency
        class(data) <- c("zooreg", "zoo")
        return(data)
    } else {
        return(zoo(data, order.by, frequency))
    }
}

rev.zooreg <- function(x) { z <- as.zooreg(rev(as.zoo(x))); frequency(z) <- frequency(x); z }
