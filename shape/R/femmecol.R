
##==============================================================================
## femmecol: red-yellow-blue colors
##==============================================================================

femmecol <- function (n=100) {

  ## red-green-blue colors on scale of 0 to 1
  rgb.col <- matrix(nrow=6,ncol=3,byrow=TRUE,
             data=c(0,0,143,
                    0,0,255,
                    0,255,255,
                    255,255,0,
                    255,0,0,
                    128,0,0))

  x.from <- c(0.0, seq(0.125,1,by=0.25), 1)  # scale from 0-1

  return(intpalette (rgb.col, n, x.from = x.from) )

}
