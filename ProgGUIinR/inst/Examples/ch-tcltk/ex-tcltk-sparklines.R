### R code from vignette source 'ex-tcltk-sparklines.Rnw'

###################################################
### code chunk number 1: ex-tcltk-sparklines.Rnw:23-30
###################################################
## sparklines
library(tcltk)
library(tseries)
window <- tktoplevel()
tkwm.title(window, "Sparklines example")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 2: ex-tcltk-sparklines.Rnw:35-40
###################################################
mL <- function(label) { # save some typing
  if(is.numeric(label))
    label <- sprintf("%.2f", label)
  ttklabel(frame, text = label, justify = "right") 
}


###################################################
### code chunk number 3: makeHeaderRule
###################################################
tkgrid(mL(""), mL("2000-01-01"), mL("-- until --"), 
       mL("today"), mL("low"), mL("high"))
tkgrid(ttkseparator(frame), row=1, column = 1, columnspan = 5, 
       sticky = "we")


###################################################
### code chunk number 4: ex-tcltk-sparklines.Rnw:56-86
###################################################
add_sparkline <- function(label, symbol = "MSFT") {
  width <- 100; height = 15               # fix width, height
  y <- get.hist.quote(instrument=symbol, start = "2000-01-01",
                      quote = "C", provider = "yahoo", 
                      retclass = "zoo")$Close
  min <- min(y); max <- max(y)
  ##
  start <- y[1]; end <- tail(y,n = 1)
  rng <- range(y)
  ##
  spark_line_canvas <- tkcanvas(frame,
                                width=width, height = height)
  x <- 0:(length(y)-1) * width/length(y)
  if(diff(rng) !=  0) {
    y1 <- (y - rng[1])/diff(rng) * height
    y1 <- height - y1   # adjust to canvas coordinates
  } else {
    y1 <- height/2 + 0 * y
  }
  ## make line with: pathName create line x1 y1... xn yn 
  l <- list(spark_line_canvas, "create","line")
  sapply(seq_along(x), function(i) {
    l[[2*i + 2]] <<- x[i]
    l[[2*i + 3]] <<- y1[i]
  })
  do.call("tcl", l)

  tkgrid(mL(label),mL(start), spark_line_canvas, 
         mL(end), mL(min), mL(max), pady = 2, sticky = "e")
}


###################################################
### code chunk number 5: ex-tcltk-sparklines.Rnw:90-93 (eval = FALSE)
###################################################
## add_sparkline("Microsoft", "MSFT")
## add_sparkline("General Electric", "GE")
## add_sparkline("Starbucks", "SBUX")


