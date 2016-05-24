## panel.barchartt is needed for S-Plus and ignored by R
## panel.barchart adds the ref= argument to the trellis function panel.barchart

panel.barchartt <-
function(x, y, col = trellis.par.get("bar.fill")$col, border = 1, ...,
         transpose = FALSE,
         ref=if (transpose) par("usr")[3] else par("usr")[1])
{
	ok <- !is.na(x) & !is.na(y)
	x <- x[ok]
	y <- y[ok]
	n <- length(y)
	NAs <- rep(NA, n)
        if (transpose) {
          ymin <- rep(ref, n)
          x1 <- x - 1/3
          x2 <- x + 1/3
          polygon(c(rbind(x1, x1, x2, x2, NAs)), c(rbind(ymin, y, y, ymin, NAs)),
                  col = col, border = as.numeric(border), ...)
          panel.abline(h=ref)
        }
        else {
          xmin <- rep(ref, n)
          y1 <- y - 1/3
          y2 <- y + 1/3
          polygon(c(rbind(xmin, x, x, xmin, NAs)), c(rbind(y1, y1, y2, y2, NAs)),
                  col = col, border = as.numeric(border), ...)
          panel.abline(v=ref)
        }
}

panel.barchart <-
if.R(r=function(...){lattice::panel.barchart(...)},
     s={
       function(..., vref=par("usr")[1]) {
         get("panel.barchart", where="trellis")(...)
         panel.abline(v=ref)
       }
     })

## panel.barchart <-
## function(x, y, col = trellis.par.get("bar.fill")$col, border = 1, ...,
##          ref=par("usr")[1])
## {
## 	ok <- !is.na(x) & !is.na(y)
## 	x <- x[ok]
## 	y <- y[ok]
## 	n <- length(y)
## 	NAs <- rep(NA, n)
##         xmin <- rep(ref, n)
##         y1 <- y - 1/3
##         y2 <- y + 1/3
##         polygon(c(rbind(xmin, x, x, xmin, NAs)), c(rbind(y1, y1, y2, y2, NAs)),
##                 col = col, border = as.numeric(border), ...)
##         panel.abline(v=ref)  ## additional line
## }
