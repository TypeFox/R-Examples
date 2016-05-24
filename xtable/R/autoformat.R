autoformat <- function(xtab, zap = getOption("digits")) {
  align(xtab) <- xalign(xtab)
  digits(xtab) <- xdigits(xtab, zap = zap)
  display(xtab) <- xdisplay(xtab)
  return(xtab)
}

xalign <- function(x, pad = TRUE) {
  lr <- function(v) if(is.numeric(v)) "r" else "l"

  is.2d <- length(dim(x)) == 2
  alignment <- if(is.2d) sapply(as.data.frame(x), lr) else lr(x)
  output <- if(is.2d && pad) c("l", alignment) else alignment

  return(output)
}

xdigits <- function(x, pad = TRUE, zap = getOption("digits")) {
  dig <- function(v) {
    if(is.numeric(v)) {
      v <- na.omit(v)
      v <- zapsmall(abs(v - floor(v)), zap)
      dec <- if(any(v > 0)) max(nchar(v) - 2L) else 0L
    } else {
      dec <- 0L
    }
    return(dec)
  }

  is.2d <- length(dim(x)) == 2
  decimals <- if(is.2d) sapply(as.data.frame(x), dig) else dig(x)
  output <- if(is.2d && pad) c(0L, decimals) else decimals

  return(output)
}

xdisplay <- function(x, pad = TRUE) {
  type <- function(v) {
    if(is.numeric(v)) {
      tp <- if(xdigits(v) == 0) "d" else "f"
    } else {
      tp <- "s"
    }
    return(tp)
  }

  is.2d <- length(dim(x)) == 2
  disp <- if(is.2d) sapply(as.data.frame(x), type) else type(x)
  output <- if(is.2d && pad) c("s", disp) else disp

  return(output)
}

