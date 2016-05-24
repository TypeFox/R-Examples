library(gWidgets2)


about <- "
This example shows the `ggraphics` component that allows one to embed a graphics device within a GUI.
The `gWidgets2tcltk` package does not support the `ggraphics` component, but one can embed a `tkrplot`
object for a similar, though different, effect.

The `ggraphics` component aims to support three mouse interactions:

* addHandlerClick to listen for user clicks (like `locator`)
* addHandlerSelectionChanged (also addHandlerChanged) to listen for the end of rubber band selection
* addHandlerMouseMotion. 

This example shows how the first two can be used to mimic 'brushing', where one graph interacts with
another
"
## Poor mans brushing
if(gtoolkit() == "tcltk")
  stop("There is no ggraphics support in tcltk. One can use `tkrplot` within a GUI though")


w <- gwindow("brushing example", visible=FALSE)
sb <- gstatusbar(cont=w)
g <- gvbox(cont=w)

pg <- gpanedgroup(cont=g, label="click, rubber band", expand=TRUE, fill=TRUE)
dev1 <- ggraphics(cont=pg)
dev2 <- ggraphics(cont=pg)

bg <- ggroup(cont=g)
addSpring(bg)
gbutton("About", cont=bg, handler=function(...) {
  w1 <- gwindow("About the example", visible=FALSE)
  g <- gvbox(cont=w1); g$set_borderwidth(10)
  glabel(about, cont=g)
  gseparator(cont=g)
  g1 <- ggroup(cont=g); addSpring(g1)
  gbutton("dismiss", cont=g1, handler=function(...) dispose(w1))
  visible(w1) <- TRUE
})

size(w) <- c(800, 300)
visible(w) <- TRUE

## Many hard coded pieces below.
x <- mtcars$mpg
y <- mtcars$wt

make_data <- function() {

  n <- 20
  bins <- seq(min(x)-1, max(x)+1, length=n+1)
  d <- cut(x, bins)
  cnts <- sapply(1:n, function(i) sum(x >= bins[i] & x < bins[i+1]))
  list(bins=bins, cnts=cnts)
}

plot_hist <- function(ind=rep(FALSE, length(x))) {
  ## make a histogram coloring some bins in

  l <- make_data()
  bins <- l$bins; cnts <- l$cnts
  n <- length(bins) - 1
  color_me <- sapply(1:n, function(i)
                     any(sapply(x[ind], function(j) j >= bins[i] & j < bins[i+1])))
  
  color <- c("gray", "red")

  visible(dev1) <- TRUE

  plot.new()
  plot.window(xlim=range(bins), ylim=range(cnts))
  mapply(rect, bins[-length(bins)], 0, bins[-1], cnts, col=color[1 + color_me])
  title("Histogram of x")
}


plot_scatter <- function(ind=rep(FALSE, length(x))) {

  visible(dev2) <- TRUE

  plot.new()
  plot.window(xlim=range(x), ylim=range(y))

  plot(y ~ x)
  points(x[ind], y[ind], col="red", pch=16, cex=2)
}


plot_hist()
plot_scatter()


## Now make interactive

## click on histogram, update scatter
addHandlerClicked(dev1, handler=function(h,...) {
  xc <- h$x
  ## which values are in x
  l <- make_data()
  bins <- l$bins; cnts <- l$cnts
  
  n <- length(bins) - 1
  which_bin <- which(sapply(1:n, function(i) bins[i] <= xc & xc < bins[i+1]))
  ind <- sapply(x,function(j) bins[which_bin] <= j & j < bins[which_bin + 1])
  plot_hist(ind)
  plot_scatter(ind)
})

addHandlerChanged(dev2, handler=function(h,...) {

  x_ind <- h$x[1] <= x & x < h$x[2]
  y_ind <- h$y[1] <= y & y < h$y[2]
  
  

  plot_hist(x_ind)
  plot_scatter(x_ind & y_ind)
})


## show which bin
addHandlerMouseMotion(dev1, handler=function(h,...) {
  xc <- h$x
  ## which values are in x
  l <- make_data()
  bins <- l$bins; cnts <- l$cnts
  
  n <- length(bins) - 1
  which_bin <- which(sapply(1:n, function(i) bins[i] <= xc & xc < bins[i+1]))
  svalue(sb) <- sprintf("In bin %s", which_bin)
})
