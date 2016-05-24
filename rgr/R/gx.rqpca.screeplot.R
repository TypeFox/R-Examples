gx.rqpca.screeplot <-
function (save, main = "", ...) 
{
    frame()
    oldpar <- par()
    on.exit(par(oldpar))
    if (main == "") 
        banner <- paste("Scree plot for", deparse(substitute(save)), 
            "\ndata source:", save$input)
    else banner <- main
    m <- length(save$eigenvalues)
    eigcum <- numeric(m)
    eigcum <- cumsum(save$econtrib)
    plot(1:m, save$eigenvalues, xlim = c(1, m), xlab = "Ordered Eigenvalues", 
        ylab = "Eigenvalue", main = banner, cex.main = 0.9, ...)
    lines(1:m, save$eigenvalues, lty = 1, ...)
    par(new = TRUE)
    plot(1:m, eigcum, axes = FALSE, xlim = c(1, m), ylim = c(1, 
        100), xlab = " ", ylab = " ", ...)
    lines(1:m, eigcum, lty = 3, ...)
    axis(side = 4)
    mtext("\nCumulative %", side = 4, line = 2, outer = FALSE)
    invisible()
}
