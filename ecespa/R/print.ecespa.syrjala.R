`print.ecespa.syrjala` <-
function(x, ...){
    cat("Syrjala test for the difference between the spatial distributions of \n",
                            x$datanames[1], "and ", x$datanames[2], ", based on", x$nsim,"simulations\n\n")
    cat("   psi:     ",x$psi.obs,"\n")
    if(mean(x$psi.sim) <= x$psi.obs)
       cat("   p-value: ", sum(x$psi.sim >= x$psi.obs)/(x$nsim+1), "\n\n")
    else cat("   p-value: ", sum(x$psi.sim <= x$psi.obs)/(x$nsim+1), "\n\n")
}

