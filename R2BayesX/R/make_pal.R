make_pal <- function(col, ncol = NULL, data = NULL, range = NULL, 
  breaks = NULL, swap = TRUE, symmetric = TRUE) 
{
  if(is.null(symmetric))
    symmetric <- TRUE
  if(is.null(col))
    col <- colorspace::diverge_hcl
  if(is.null(ncol) && is.null(breaks))
    ncol <- 99L
  if(is.null(ncol) && !is.null(breaks))
    ncol <- length(breaks) - 1L
  if(is.function(col))
    col <- col(ncol)    
  else 
    ncol <- length(col)
  if(swap) 
    col <- rev(col)
  if(all(is.null(data), is.null(range), is.null(breaks))) 
    stop("at least one needs to be specified")
  if(is.null(breaks)) {
    if(is.null(range)) {
      range <- range(data, na.rm = TRUE)
      if(symmetric) { 
        mar <- max(abs(range))
        range <- c(0 - mar, mar)
      }
    }
    if(diff(range) == 0)
      breaks <- seq(min(range) - 1, min(range) + 1, length.out = ncol + 1L)
    else
      breaks <- seq(range[1L], range[2L], length.out = ncol + 1L)
  } else stopifnot(length(breaks) == ncol + 1L)
  if(is.matrix(data)) {
    obs2col <- function(x) {
      hgt <- (x[-1L, -1L] + x[-1L, -(ncol(x) - 1L)] + 
        x[-(nrow(x) -1L), -1L] + x[-(nrow(x) -1L), -(ncol(x) - 1L)])/4
      c(col[1L], col, col[ncol])[cut(hgt, c(-Inf, breaks, Inf))]
      }
  } else {
    obs2col <- function(x) c(col[1L], col, col[ncol])[cut(x, c(-Inf, breaks, Inf))]
  }

  return(list(colors = col, breaks = breaks, map = obs2col))
}

