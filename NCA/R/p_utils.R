p_generate_title <- 
function (loop.data) {
  return (paste(loop.data$names[1], loop.data$names[2], sep = " - "))
}

p_pretty_name <- 
function (uglyName) {
  return( gsub("_", "-", toupper(uglyName)) )
}

p_pretty_number <- 
function (uglyNumber, default="", prec=3, useSpaces=FALSE) {
  if (is.na(uglyNumber) || uglyNumber == "NA") {
    return(default)
  } 
  
  if (is.integer(uglyNumber) && !useSpaces) {
    return(sprintf("%d", uglyNumber))
  }
  
  if (prec == "auto") {
    if (uglyNumber == 0) { 
      prec <- 3
    } else {
      prec <- max(0, 3 - floor(log10(abs(uglyNumber))))
    }
  }
  fmt <- sprintf("%%.%df%%s", prec)
  
  nSpaces <- 0
  if (useSpaces) {
    nSpaces <- ifelse(prec == 0, 4, max(0, 3-prec))
  }

  # We hate to see -0.0
  uglyNumber[abs(uglyNumber) < 0.1 ^ prec] <- 0

  return(sprintf(fmt, uglyNumber, paste(rep(" ", nSpaces), collapse='')))
}

p_pretty_mpx <- 
function (loop.data, mpx, max.value) {
  if (all(is.na(mpx))) {
    colnames(mpx) <- c(loop.data$names[1])
    return (mpx)    
  }
  
  mpxMax <- max(mpx, na.rm=TRUE)
  maxId  <- which.max(mpx==mpxMax)
  useMax <- function (x, y) {
    if (y > maxId && is.na(x)) {
      return (max.value)
    }
    return (x)
  }

  mpx <- mapply(useMax, mpx, row(mpx))
  mpx <- matrix(mpx, nrow=length(mpx))
  colnames(mpx) <- c(loop.data$names[1])
  
  return(mpx)
}

p_warn_percentage_max <-
function (bottleneck.y, loop.data, id.y) {
  if (p_bottleneck_id(bottleneck.y) == 2 && loop.data$y.low < 0) {
    message("")
    message(paste0("Warning : using bottleneck.y with Y values < 0",
                   ", results might be counterintuitive!"))
    message("")
  }
}
