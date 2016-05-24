display.ascii.o <-
function () 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    par(usr = c(0, 17, 0, 17))
    mtext("Windows 3.1 Latin 1 Octal Table\nUse in R text strings, e.g., \\265 for mu", 
        side = 1, line = 1.1, cex = 1.2)
    on <- structure(0:256, class = "octmode")
    for (y in 1:16) for (x in 1:16) {
        n <- 16 * (y - 1) + x - 1
        points(x + 0.1, y, pch = n, cex = 0.7)
        text(x - 0.3, y + 0.5, paste(on[n + 1]), cex = 0.6, adj = 0)
    }
    lines(c(0.1, 0.1, 16.9, 16.9, 0.1), c(0.5, 17, 17, 0.5, 0.5))
    invisible()
}
