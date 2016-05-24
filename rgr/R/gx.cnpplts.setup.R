gx.cnpplts.setup <-
function (display = FALSE) 
{
    pchs <- c(0, 5, 2, 6, 1, 10, 12, 13, 14)
    symcols <- c(1, 2, 4, 3, 5, 1, 6, 4, 3)
    cex <- 0.8
    cexp <- 0.9
    if (display) {
        frame()
        oldpar <- par()
        on.exit(par(oldpar))
        par(usr = c(0, 3, 0, 3))
        mtext("Symbology for function gx.cnpplts", side = 3, 
            line = 2, cex = 2)
        for (i in 1:9) {
            x <- ((i - 1)%%3) + 0.5
            y <- 2.5 - ((i - 1)%/%3)
            points(x, y, pch = pchs[i], col = symcols[i], cex = 7)
            text(x, y + 0.35, paste("subset =", i), adj = 0.5)
            text(x, y - 0.35, paste("pchs =", pchs[i], "\nsymcols =", 
                symcols[i]), adj = 0.5)
        }
        abline(h = 0:3, lty = 1)
        abline(v = 0:3, lty = 1)
    }
    invisible(list(pchs = pchs, symcols = symcols, cex = cex, 
        cexp = cexp))
}
