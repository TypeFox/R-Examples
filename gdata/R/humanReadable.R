humanReadable <- function(x,
                          units="auto",
                          standard=c("IEC", "SI", "Unix"),
                          digits=1,
                          width=NULL,
                          sep=" ",
                          justify = c("right", "left")
                          )
{
  ## --- Setup ---

  suffix.SI   <- c("B",  "kB",  "MB",  "GB",  "TB",  "PB",  "EB",  "ZB",  "YB")
  suffix.IEC  <- c("B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB")
  suffix.Unix <- c("B" ,  "K",   "M",   "G",   "T",   "P",   "E",   "Z",   "Y")

  standard <- match.arg(standard)
  if(length(justify)==1) justfy <- c(justify, justify)

  ## --- Functions ---

  .applyHuman <- function(x, base, suffix, digits, width, sep)
  {
    ## Which suffix should we use?
    n <- length(suffix)
    i <- pmax(pmin(floor(log(x, base)), n-1),0)
    if(!is.finite(i)) i <- 0
    x <- x / base^i
    ## Formatting
    if(is.null(width))
        ## the same formatting for all
        x <- format(round(x=x, digits=digits), nsmall=digits)
    else
        {
            ## similar to ls, du, and df
            lenX <- nchar(x)
            if(lenX > width) {
                digits <- pmax( width - nchar(round(x)) - 1, 0)
            }
            if(i == 0) digits <- 0
            x <- round(x, digits=digits)
        }
    c(x, suffix[i+1])
  }

  ## -- Work

  if(any(x < 0)) stop("'x' must be positive")
  if(standard == "SI")
      {
          suffix <- suffix.SI
          base <- 10^3
      }
  else if (standard=="IEC")
      {
          suffix <- suffix.IEC
          base <- 2^10
      }
  else # (standard=="Unix)
      {
          suffix <- suffix.Unix
          base <- 2^10
      }

  if(!missing(units) && units=="bytes")
      {
          retval <- rbind(x, "bytes")
      }
  else if(!missing(units) && units!="auto")
      {
          units <- suffix[match( toupper(units), toupper(suffix) )]
          power <- match(units, suffix ) -1
          X <- x/(base^power)
          X <- format.default(x=X, digits=digits, nsmall=digits)
          retval <- rbind(X, rep(units, length(X)))
      }
  else
      retval <- sapply(X=x, FUN=".applyHuman", base=base, suffix=suffix,
                       digits=digits, width=width, sep=sep)

  if(all(justify == "none"))
      paste(trim(retval[1,]), trim(retval[2,]), sep=sep)
  else
      paste(format(trim(retval[1,]), justify=justify[1]),
            format(trim(retval[2,]), justify=justify[2]),
            sep=sep)

}

