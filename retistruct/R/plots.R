##
## Utilities for plotting functions
## 

##' @title Generate colours for strain plots
##' @param x Vector of values of log strain
##' @return Vector of colours corresponding to strains
##' @author David Sterratt
strain.colours <- function(x) {
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  col <- x/log(0.75)*35 + 35
  col[col<1] <- 1
  col[col>70] <- 70
  return(col)
}

##' Place text at bottom right of \code{\link{projection}}
##'
##' @title Put text on the polar plot
##' @param text Test to place
##' @author David Sterratt
##' @export
polartext <- function(text) {
  mtext(text, 1, adj=1, line=-1)
}

##' @title Ancillary function to place labels
##' @param panlabel Label text
##' @param line Line on which to appear
##' @author David Sterratt
##' @export
panlabel <- function(panlabel, line=-0.7) {
  mtext(panlabel, side=3, adj=-par("plt")[1]/(par("plt")[2]-par("plt")[1]),
        line=line, font=2, cex=1)
}
