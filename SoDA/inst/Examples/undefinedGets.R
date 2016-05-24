`undefined<-` <- function(x, codes = numeric(), value) {
      if(length(codes) > 0)
        x[ x %in% codes] <- NA
      x[is.na(x)] <- value
      x
  }
