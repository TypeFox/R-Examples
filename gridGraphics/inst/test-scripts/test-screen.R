
library(gridGraphics)

screen1 <- function() {
    par(bg = "white")           # default is likely to be transparent
    split.screen(c(2, 1))       # split display into two screens
    split.screen(c(1, 3), screen = 2) # now split the bottom half into 3
    screen(1) # prepare screen 1 for output
    plot(10:1)
    screen(4) # prepare screen 4 for output
    plot(10:1)
    close.screen(all = TRUE)    # exit split-screen mode
}

screen2 <- function() {
    split.screen(c(2, 1))       # split display into two screens
    split.screen(c(1, 2), 2)    # split bottom half in two
    plot(1:10)                  # screen 3 is active, draw plot
    erase.screen()              # forgot label, erase and redraw
    plot(1:10, ylab = "ylab 3")
    screen(1)                   # prepare screen 1 for output
    plot(1:10)
    screen(4)                   # prepare screen 4 for output
    plot(1:10, ylab = "ylab 4")
    screen(1, FALSE)            # return to screen 1, but do not clear
    plot(10:1, axes = FALSE, lty = 2, ylab = "")  # overlay second plot
    axis(4)                     # add tic marks to right-hand axis
    title("Plot 1")
    close.screen(all = TRUE)    # exit split-screen mode
}

plotdiff(expression(screen1()), "screen-1", width=10)
plotdiff(expression(screen2()), "screen-2")

plotdiffResult()
