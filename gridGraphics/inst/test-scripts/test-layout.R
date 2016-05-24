
library(gridGraphics)

layout1 <- function() {
    ## divide the device into two rows and two columns
    ## allocate figure 1 all of row 1
    ## allocate figure 2 the intersection of column 2 and row 2
    layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE))
    ## show the regions that have been allocated to each plot
    layout.show(2)
}

layout2 <- function() {
    ## divide device into two rows and two columns
    ## allocate figure 1 and figure 2 as above
    ## respect relations between widths and heights
    nf <- layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE), respect = TRUE)
    layout.show(nf)
}

layout3 <- function() {
    ## create single figure which is 5cm square
    nf <- layout(matrix(1), widths = lcm(5), heights = lcm(5))
    layout.show(nf)
}

layout4 <- function() {
    set.seed(1)
    x <- pmin(3, pmax(-3, stats::rnorm(50)))
    y <- pmin(3, pmax(-3, stats::rnorm(50)))
    xhist <- hist(x, breaks = seq(-3,3,0.5), plot = FALSE)
    yhist <- hist(y, breaks = seq(-3,3,0.5), plot = FALSE)
    top <- max(c(xhist$counts, yhist$counts))
    xrange <- c(-3, 3)
    yrange <- c(-3, 3)
    nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(3,1), c(1,3), TRUE)
    # layout.show(nf)
    par(mar = c(3,3,1,1))
    plot(x, y, xlim = xrange, ylim = yrange, xlab = "", ylab = "")
    par(mar = c(0,3,1,1))
    barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)
    par(mar = c(3,0,1,1))
    barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0,
            horiz = TRUE)
}

plotdiff(expression(layout1()), "layout-1")
plotdiff(expression(layout2()), "layout-2")
plotdiff(expression(layout3()), "layout-3")
plotdiff(expression(layout4()), "layout-4")

plotdiffResult()
