plot.syrjala.test <-
function (x,coline = 1, ...) 
{
    hist(x$cvm.sim, col = "grey", border = "grey", main = "", 
        xlab = expression(psi), ...)
    title(main = paste("Cramer-von Misses test\nfor the difference between the spatial distributions of\n", 
        x$datanames[1], "and", x$datanames[2]))
    abline(v = x$cvm.obs, lwd = 3,col=coline)
}

