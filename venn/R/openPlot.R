`openPlot` <-
function(size = 15) {
    if (dev.cur() == 1) {
        dev.new(width = (size + 1)/2.54, height = (size + 1)/2.54)
    }
    
    par(new = FALSE, xpd = TRUE, mai = c(0.05, 0.05, 0.05, 0.05))
    plot(0:1000, type = "n", axes = FALSE, asp = 1)
}

