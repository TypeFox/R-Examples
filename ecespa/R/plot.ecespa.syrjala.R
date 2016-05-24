`plot.ecespa.syrjala` <-
function(x, ...){
    hist(x$psi.sim, col="grey", border="grey", main="", xlab=expression(psi),...)
    title(main=paste("Syrjala test\n for the difference between the spatial distributions of \n",
                            x$datanames[1], " and ", x$datanames[2]))
    abline(v=x$psi.obs, lwd=3)
    print(x)
}

