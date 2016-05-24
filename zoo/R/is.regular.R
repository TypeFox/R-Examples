is.regular <- function(x, strict = FALSE) {
  UseMethod("is.regular")
}

is.regular.zoo <- function(x, strict = FALSE)
{
  delta <- suppressWarnings(try(diff(as.numeric(index(x))), silent = TRUE))
  if(class(delta) == "try-error" || any(is.na(delta))) FALSE
  else if(length(delta) < 1) FALSE
  else if(strict) identical(all.equal(delta, rep.int(delta[1], length(delta))), TRUE)
  else {
    delta <- unique(delta)
    rval <- identical(all.equal(delta/min(delta), round(delta/min(delta))), TRUE)
    if(!rval && identical(all.equal(delta, round(delta)), TRUE)) rval <- TRUE
    rval
  }
}

is.regular.ts <- function(x, strict = FALSE) TRUE

is.regular.zooreg <- function(x, strict = FALSE)
{
  if(strict) is.regular.zoo(x, strict = TRUE) else TRUE
}

is.regular.default <- function(x, strict = FALSE) {
  is.regular(as.zoo(x), strict = strict)
}

frequency.zooreg <- function(x, ...) 
{
  attr(x, "frequency")
}

frequency.zoo <- function(x, ...)
{
  ## check whether frequency is available
  freq <- attr(x, "frequency")
  if(!is.null(freq) || length(index(x)) < 2) return(freq)

  ## check regularity
  delta <- suppressWarnings(try(diff(as.numeric(index(x))), silent = TRUE))
  reg <- if(class(delta) == "try-error" || any(is.na(delta))) FALSE
  else {
    delta <- unique(delta)
    rval <- identical(all.equal(delta/min(delta), round(delta/min(delta))), TRUE)
    if(rval) freq <- 1/min(delta)
    else if(identical(all.equal(delta, round(delta)), TRUE)) {
      ## special case: integer indexes
      ## get frequency as greatest common divisor (of differences)
      gcd <- function(x) {	
        gcd0 <- function(a, b) ifelse(b==0 | a==b, a, gcd0(b, a %% b))
        if(length(x) < 2) x <- c(x, as.integer(0))
        if(length(x) < 3) {
          return(gcd0(x[1], x[2]))
        } else {
          x <- sapply(1:(length(x) - 1), function(i) gcd0(x[i], x[i+1]))
          gcd(x)
        }
      }
      freq <- 1/gcd(delta)
      rval <- TRUE
    }
    rval
  }
  if(!reg) return(NULL)
  if(freq > 1 && identical(all.equal(freq, round(freq)), TRUE)) freq <- round(freq)
  return(freq)
}

"frequency<-" <- function(x, value)
  UseMethod("frequency<-")
  
"frequency<-.zoo" <- function(x, value) {
  delta <- suppressWarnings(try(diff(as.numeric(index(x))), silent = TRUE))
  freqOK <- if(class(delta) == "try-error" || any(is.na(delta))) FALSE
    else if(length(delta) < 1) TRUE
    else identical(all.equal(delta*value, round(delta*value)), TRUE)
  stopifnot(freqOK)
  if(value > 1 && identical(all.equal(value, round(value)), TRUE)) value <- round(value)
  attr(x, "frequency") <- value
  class(x) <- c("zooreg", "zoo")
  return(x)
}

"frequency<-.zooreg" <- function(x, value) {
  delta <- diff(as.numeric(index(x)))
  stopifnot(identical(all.equal(delta*value, round(delta*value)), TRUE))
  attr(x, "frequency") <- value
  return(x)
}

deltat.zoo <- function(x, ...)
{
  rval <- frequency.zoo(x, ...)
  if(is.null(rval)) NULL else 1/rval
}

deltat.zooreg <- function(x, ...)
{
  1/frequency.zooreg(x, ...)
}

cycle.zooreg <- function(x, ...)
{
  freq <- frequency(x)
  ix <- as.numeric(index(x))
  d <- diff(ix)
  if(!identical(all.equal(freq*d, round(freq*d)), TRUE))
    stop(paste(sQuote("cycle"), "not available for", sQuote("x")))  
  return(zoo(round((ix - floor(ix)) * freq) + 1, order.by = index(x), freq))
}

cycle.zoo <- function(x, ...)
{
  if(is.regular(x)) cycle.zooreg(x)
    else stop(sQuote("x"), "is not regular")
}
