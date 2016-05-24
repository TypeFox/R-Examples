
##==============================================================================
## greycol: black-white colors
##==============================================================================

greycol <- function (n=100, interval = c(0.0,0.7)) {

  return(shadepalette (n=n, inicol="white", endcol="black", interval=interval))
}

# alias for greycol...
graycol <- function (n=100,                 # number of colors
                     interval = c(0.0,0.7)) # interval *to* where to interpolate
  return(shadepalette(n=n,inicol="white",endcol="black",interval=interval))
