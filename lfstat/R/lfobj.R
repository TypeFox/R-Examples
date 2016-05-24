if(getRversion() >= "2.15.1"){
  utils::globalVariables(c("day", "month", "year", "flow", "tp", "baseflow",
                           ".warned"), add = TRUE)
}

#Different methods to create a lfobj:
#Data.frame with named columns
#ts + start date (does vector work?)

createlfobj <- function(x, ...){
  UseMethod("createlfobj")
}


createlfobj.lfobj <- function(x, hyearstart = NULL, baseflow = NULL,
                              meta = NULL, ...){
  if(is.null(baseflow)){
    baseflow <- "baseflow" %in% names(x)
  }

  if(is.null(meta)){
    meta <- attr(x, "lfobj")
  }

  dat <- createlfobj.data.frame(x = x, hyearstart = hyearstart,
                                baseflow = baseflow,
                                meta = meta, ...)
  return(dat)
}


# Create a lfobj from a vector of daily flow data and the startdate
createlfobj.ts <- function(x, startdate, dateformat = "%d/%m/%Y", ...){

  start <- as.Date(startdate, dateformat)
  time <- seq(from = start, along.with = x, by = "days")
  df <- data.frame(strsplit_date(time), flow = as.vector(x))

  dat <- createlfobj(x = df, ...)
  return(dat)
}


#Create a lfobj from a data frame with cols named "flow", "day", "month", "year"
createlfobj.data.frame <- function(x, hyearstart = NULL, baseflow = TRUE,
                                   meta = list(), ...){

  cols <- c("day", "month", "year", "flow")
  if(!all(cols %in% names(x))) {
    stop("Your data frame must contain colums named",
         paste(shQuote(cols), collapse = ", "),
         "! Please look at the help files for more information.")
  }


  if(!(is.null(hyearstart) || hyearstart %in% 1:12)){
    stop("if set, hyearstart must be an integer between 1 and 12")
  }


  # try to guess from column hyear
  if((is.null(hyearstart) || (!hyearstart %in% 1:12))){
    hyearstart <- hyear_start(x)
  }



  meta <- as.list(meta)
  meta[["hyearstart"]] <- hyearstart
  x <- as.data.frame(x)

  dat <- x[, cols]
  time <- with(x, as.Date(paste(year, month, day, sep = "-")))

  # hydrological year is kept as numeric for backwards compatibility
  dat$hyear <- as.numeric(as.character(water_year(time, origin = hyearstart)))


  fullseq <- seq(from = min(time), to = max(time), by = "day")
  missing <- fullseq[!fullseq %in% time]
  if(length(missing)) {
    warning("Irregular time series provided. Missing obervations were padded with NAs.")
    gaps <- data.frame(strsplit_date(missing), flow = NA)
    dat <- rbind(dat, gaps)
  }

  # reorder if nescessary
  if(is.unsorted(time) || length(missing)) dat <- dat[order(c(time, missing)), ]
  rownames(dat) <- NULL

  if(baseflow) dat$baseflow <- baseflow(dat$flow, ...)

  # Meta-Information
  attr(dat, "lfobj") <- meta

  class(dat) <- c("lfobj", "data.frame")
  return(dat)
}


# hack to make attributes sticky
# otherwise subsetting would loose attributes
"[.lfobj" <- function (x, i, j, drop = T) {

  y <- "[.data.frame"(x, i, j, drop)
  attr(y, "lfobj") <- attr(x, "lfobj")

  return(y)
}





lfcheck <- function(lfobj){
  if(!inherits(lfobj,"lfobj")){
    stop("This functions is designed for objects of the class 'lfobj'. ",
         "Please use 'createlfobj()' or see '?createlfobj' for more information")
  }
}


